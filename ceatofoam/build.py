"""From an input dictionary to a design and an OpenFOAM case.

Input (YAML or a dict); quantities are numbers in SI or strings with units::

    propellants: {fuel: RP-1, oxidizer: O2(L)}     # optional temperatures:
                                                   # fuel_temperature, ...
    chamber:     {pressure: 500 psi, mixture_ratio: 2.4}
    design:      {thrust: 500 lbf, altitude: 0 ft, sizing: frozen}
    geometry:    {characteristic_length: 40 in,
                  convergent_half_angle: 30 deg, divergent_half_angle: 15 deg}
    mesh:        {radial_cells: 40}
    solver:      {model: kOmegaSST, iterations: 20000}

Instead of ``propellants`` a fixed-composition gas can be given (cold-flow
validation cases), with the chamber temperature::

    gas:     {composition: {N2: 0.78084, O2: 0.20946, AR: 0.00934},
              mechanism: air.yaml}
    chamber: {pressure: ..., temperature: ...}

``design.sizing`` is frozen (default), equilibrium or constant-gamma
(``design.gamma``: a number, ``gammas`` for the CEA chamber gamma_s or
``frozen`` for the frozen chamber gamma; ``design.molecular_weight``:
``chamber`` (default) or a number). The ambient state is either
``design.altitude`` or ``design.ambient_pressure`` and
``design.ambient_temperature``.

Geometry: either sized from the design (``characteristic_length``,
``contraction_ratio`` [auto: Humble et al. correlation]) or explicit
(``throat_radius``, ``area_ratio``, ``contraction_ratio``,
``chamber_length``). Arc radii are multiples of the throat radius:
``throat_upstream_radius`` and ``throat_downstream_radius`` (default 1.5,
as in ceaToFoam 1.x) and ``chamber_fillet_radius`` (default 1.0).
"""

from dataclasses import asdict
import datetime
import json
import math
import os
import platform
import subprocess

import numpy as np

from . import __version__, atmosphere, cea, geometry, mesh, sizing, thermo
from .foamcase import CaseConditions, SolverSettings, write_case
from .units import to_si


def load(path):
    import yaml

    with open(path) as f:
        return yaml.safe_load(f)


def _ambient(design):
    if "ambient_pressure" in design:
        return (
            to_si(design.get("ambient_temperature", 288.15), "temperature"),
            to_si(design["ambient_pressure"], "pressure"),
        )
    return atmosphere.ambient(to_si(design.get("altitude", 0.0), "length"))


def _subsonic_mach(area_ratio, gamma):
    """Subsonic Mach number at A/A* = area_ratio (constant gamma)."""
    from scipy.optimize import brentq

    return brentq(
        lambda m: sizing.area_ratio_from_mach(m, gamma) - area_ratio,
        1e-6, 1.0 - 1e-9,
    )


def make_design(cfg):
    """Return a dict with the chamber data, mixture, design and contour."""
    chamber = cfg["chamber"]
    design_cfg = cfg.get("design", {})
    geo = cfg.get("geometry", {})
    pc = to_si(chamber["pressure"], "pressure")
    t_amb, p_amb = _ambient(design_cfg)
    model = design_cfg.get("sizing", "frozen")

    cea_points = {}
    if "propellants" in cfg:
        prop = cfg["propellants"]
        of = float(chamber["mixture_ratio"])
        kwargs = dict(
            t_fuel=prop.get("fuel_temperature"),
            t_oxidizer=prop.get("oxidizer_temperature"),
        )
        if kwargs["t_fuel"] is not None:
            kwargs["t_fuel"] = to_si(kwargs["t_fuel"], "temperature")
        if kwargs["t_oxidizer"] is not None:
            kwargs["t_oxidizer"] = to_si(kwargs["t_oxidizer"], "temperature")
        for analysis in ("equilibrium", "frozen"):
            cea_points[analysis] = cea.run_rocket(
                prop["fuel"], prop["oxidizer"], pc, of, pc/p_amb, analysis,
                **kwargs,
            )
        tc = cea_points["frozen"].tc
        mixture = thermo.FrozenMixture(
            cea_points["frozen"].composition,
            cfg.get("mechanism", "gri30_highT.yaml"),
        )
    elif "gas" in cfg:
        gas = cfg["gas"]
        tc = to_si(chamber["temperature"], "temperature")
        mixture = thermo.FrozenMixture(
            gas["composition"], gas.get("mechanism", "gri30_highT.yaml")
        )
    else:
        raise ValueError("input needs 'propellants' or 'gas'")

    explicit = "throat_radius" in geo
    thrust = (
        to_si(design_cfg["thrust"], "force") if "thrust" in design_cfg
        else None
    )

    if model == "frozen":
        target = sizing.size_frozen(
            thrust if thrust else 1.0, pc, tc, p_amb, mixture
        )
    elif model == "equilibrium":
        if "equilibrium" not in cea_points:
            raise ValueError("equilibrium sizing needs propellants")
        target = sizing.size_equilibrium(
            thrust if thrust else 1.0, p_amb, cea_points["equilibrium"]
        )
    elif model == "constant-gamma":
        g = design_cfg.get("gamma", "gammas")
        if g == "gammas":
            if "equilibrium" not in cea_points:
                raise ValueError("gamma: gammas needs propellants")
            gamma = cea_points["equilibrium"].gammas_c
        elif g == "frozen":
            gamma = mixture.gamma(tc)
        else:
            gamma = float(g)
        mw = design_cfg.get("molecular_weight", "chamber")
        mw = mixture.mw if mw == "chamber" else float(mw)
        target = sizing.size_constant_gamma(
            thrust if thrust else 1.0, pc, tc, p_amb, gamma, mw
        )
    else:
        raise ValueError(f"unknown sizing model {model!r}")

    if explicit:
        rt = to_si(geo["throat_radius"], "length")
        eps = float(geo["area_ratio"])
        ec = float(geo["contraction_ratio"])
        # The reference (ideal) mass flow of this throat with frozen flow.
        throat = sizing.frozen_throat(mixture, tc, pc)
        mdot_ideal = throat.mass_flux*math.pi*rt**2
    else:
        if thrust is None:
            raise ValueError("sized geometry needs design.thrust")
        rt = target.rt
        eps = target.area_ratio
        ec = geo.get("contraction_ratio", "auto")
        ec = (
            geometry.contraction_ratio_correlation(2*rt) if ec == "auto"
            else float(ec)
        )
        mdot_ideal = target.mdot

    rc = rt*math.sqrt(ec)
    re = rt*math.sqrt(eps)
    at = math.pi*rt**2
    if explicit and "chamber_length" in geo:
        # Explicit cylinder length: chamber volume follows from it.
        l_cyl = to_si(geo["chamber_length"], "length")
        probe = geometry.conical_contour(
            rt, rc, re,
            to_si(geo.get("convergent_half_angle", 30), "angle"),
            to_si(geo.get("divergent_half_angle", 15), "angle"),
            float(geo.get("chamber_fillet_radius", 1.0))*rt,
            float(geo.get("throat_upstream_radius", 1.5))*rt,
            float(geo.get("throat_downstream_radius", 1.5))*rt,
            1e9,
        )
        v_conv = probe.volume(
            probe.segments[1].start[0], probe.x_throat
        )
        volume = v_conv + math.pi*rc**2*l_cyl
    else:
        lstar = to_si(geo.get("characteristic_length", "40 in"), "length")
        volume = lstar*at

    contour = geometry.conical_contour(
        rt, rc, re,
        to_si(geo.get("convergent_half_angle", 30), "angle"),
        to_si(geo.get("divergent_half_angle", 15), "angle"),
        float(geo.get("chamber_fillet_radius", 1.0))*rt,
        float(geo.get("throat_upstream_radius", 1.5))*rt,
        float(geo.get("throat_downstream_radius", 1.5))*rt,
        volume,
    )
    return {
        "pc": pc, "tc": tc, "p_amb": p_amb, "t_amb": t_amb,
        "cea": cea_points, "mixture": mixture, "target": target,
        "explicit_geometry": explicit, "mdot_ideal": mdot_ideal,
        "contour": contour,
    }


def initial_slices(mixture, tc, pc, contour, n_sub=30, n_sup=30):
    """1-D frozen state in axial slices of the chamber: (x0, x1, p, T, u).

    Upstream of the throat the subsonic branch, downstream the supersonic
    branch, each at the area ratio of the slice centre.
    """
    throat = sizing.frozen_throat(mixture, tc, pc)
    edges_sub = np.linspace(0.0, contour.x_throat, n_sub + 1)
    edges_sup = np.linspace(contour.x_throat, contour.x_exit, n_sup + 1)
    out = []
    for edges, sup in ((edges_sub, False), (edges_sup, True)):
        for x0, x1 in zip(edges[:-1], edges[1:]):
            xm = 0.5*(x0 + x1)
            eps = float((contour.radius(xm)/contour.rt)**2)
            st = sizing.frozen_station_at_area(
                mixture, tc, pc, eps, sup, throat
            )
            # Small overlap so that no cell centre falls between boxes.
            out.append((x0 - 1e-9, x1 + 1e-9, st.p, st.T, st.u))
    return out


def _git_revision():
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        out = subprocess.run(
            ["git", "-C", here, "describe", "--always", "--dirty",
             "--abbrev=12"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def _versions():
    v = {"ceatofoam": __version__, "source_revision": _git_revision(),
         "python": platform.python_version()}
    for mod in ("cantera", "CEA_Wrap", "numpy", "scipy"):
        try:
            v[mod] = __import__(mod).__version__
        except (ImportError, AttributeError):
            v[mod] = None
    v["openfoam"] = os.environ.get("WM_PROJECT_VERSION")
    return v


def make_case(cfg, case_dir):
    """Write the OpenFOAM case for ``cfg`` into ``case_dir``.

    Returns the provenance record, which is also written to
    ``case_dir/ceatofoam.json``.
    """
    d = make_design(cfg)
    mix, contour, target = d["mixture"], d["contour"], d["target"]
    mesh_s = mesh.MeshSettings(**cfg.get("mesh", {}))
    solver_s = SolverSettings(**cfg.get("solver", {}))

    janaf = thermo.fit_janaf(mix)
    t_lo = 0.8*min(target.te, d["t_amb"]) if "gas" in cfg else 0.8*target.te
    suth = thermo.fit_sutherland(mix, t_lo, 1.05*d["tc"])

    gamma_c = mix.gamma(d["tc"])
    m_in = _subsonic_mach(contour.contraction_ratio, gamma_c)
    t_in = d["tc"]/(1 + 0.5*(gamma_c - 1)*m_in**2)
    u_in = m_in*math.sqrt(gamma_c*mix.R*t_in)

    slices = initial_slices(mix, d["tc"], d["pc"], contour)

    os.makedirs(case_dir, exist_ok=True)
    mesh_summary = mesh.write_block_mesh_dict(case_dir, contour, mesh_s)
    conditions = CaseConditions(
        pc=d["pc"], tc=d["tc"], p_amb=d["p_amb"], t_amb=d["t_amb"],
        mw=mix.mw, gamma_c=gamma_c, gamma_amb=mix.gamma(d["t_amb"]),
        janaf=janaf, sutherland=suth, x_throat=contour.x_throat,
        x_exit=contour.x_exit, rc=contour.rc, re=contour.re,
        inlet_velocity=u_in, initial_slices=slices,
        x_end=mesh_summary["x_end"],
        wall_points=[
            (float(x), float(contour.radius(x)))
            for x in np.linspace(0.0, contour.x_exit, 201)[1:-1]
        ],
    )
    write_case(case_dir, conditions, solver_s)

    record = {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "versions": _versions(),
        "input": cfg,
        "ambient": {"p": d["p_amb"], "T": d["t_amb"]},
        "chamber": {"p": d["pc"], "T": d["tc"]},
        "cea": {k: asdict(v) for k, v in d["cea"].items()},
        "mixture": {
            "mechanism": mix.mechanism,
            "composition": mix.composition,
            "dropped": mix.dropped,
            "dropped_fraction": mix.dropped_fraction,
            "molecular_weight": mix.mw,
            "gamma_chamber": gamma_c,
            "gamma_ambient": mix.gamma(d["t_amb"]),
        },
        "fits": {
            "janaf": {
                "low": list(janaf.low), "high": list(janaf.high),
                "cp_max_rel_error": janaf.cp_max_rel_error,
                "h_max_abs_error": janaf.h_max_abs_error,
            },
            "sutherland": asdict(suth),
        },
        "design": target.as_dict(),
        "explicit_geometry": d["explicit_geometry"],
        "mdot_ideal": d["mdot_ideal"],
        "geometry": {
            "rt": contour.rt, "rc": contour.rc, "re": contour.re,
            "area_ratio": contour.area_ratio,
            "contraction_ratio": contour.contraction_ratio,
            "theta_c": contour.theta_c, "theta_d": contour.theta_d,
            "r_cf": contour.r_cf, "r_tu": contour.r_tu, "r_td": contour.r_td,
            "cylinder_length": contour.cylinder_length,
            "x_throat": contour.x_throat, "x_exit": contour.x_exit,
            "chamber_volume": contour.volume(0.0, contour.x_throat),
            "stations": contour.stations,
        },
        "inlet": {"mach": m_in, "velocity": u_in},
        "mesh": {**asdict(mesh_s), **mesh_summary},
        "solver": asdict(solver_s),
    }
    with open(os.path.join(case_dir, "ceatofoam.json"), "w") as f:
        json.dump(record, f, indent=2, default=float)
    return record
