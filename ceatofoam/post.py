"""Thrust, mass flow and exit state from a finished case.

Wedge integrals are scaled to the full annulus with the ratio of the
circular area pi r^2 to the measured patch area (the same ratio holds at
every radius of a planar-faced wedge).

Thrust is evaluated in two ways:

* exit plane: F = sum(phi U_x) + integral (p - p_amb) dA over the exit faceZone;
* wall budget: F = [sum(phi U_x) + integral (p - p_amb) dA] at the injector
  face minus the axial force of the gas on the wall.

Both are equal if axial momentum is conserved, so their difference measures
the momentum budget closure of the discretization. Steadiness is the
relative range of the exit mass flow and thrust over the last
``window`` fraction of the iterations.
"""

import json
import math
import os

import numpy as np

from .units import G0


def _read_dat(path):
    """Read an OpenFOAM function-object .dat file into (names, array)."""
    names, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                tokens = line[1:].split()
                if tokens and tokens[0] == "Time":
                    names = tokens
                continue
            line = line.replace("(", " ").replace(")", " ")
            if line.strip():
                rows.append([float(v) for v in line.split()])
    return names, np.array(rows)


def _series(case, fo, filename="surfaceFieldValue.dat"):
    base = os.path.join(case, "postProcessing", fo)
    if not os.path.isdir(base):
        raise FileNotFoundError(f"no function-object output {base}")
    # Concatenate restarts (time directories), later ones overriding.
    data = None
    for t in sorted(os.listdir(base), key=float):
        path = os.path.join(base, t, filename)
        if os.path.exists(path):
            names, arr = _read_dat(path)
            if arr.size == 0:
                continue
            data = arr if data is None else np.vstack(
                [data[data[:, 0] < arr[0, 0]], arr]
            )
    if data is None:
        raise FileNotFoundError(f"no data in {base}")
    return data


def evaluate(case, window=0.1):
    with open(os.path.join(case, "ceatofoam.json")) as f:
        rec = json.load(f)
    p_amb = rec["ambient"]["p"]
    rc, re = rec["geometry"]["rc"], rec["geometry"]["re"]

    flux_e = _series(case, "exitFlux")      # Time, area, sum(phi)
    pres_e = _series(case, "exitPressure")  # Time, areaIntegrate(p)
    mom_e = _series(case, "exitMomentum")   # Time, Ux Uy Uz
    avg_e = _series(case, "exitAverage")    # Time, p T Ma
    flux_i = _series(case, "inletFlux")
    pres_i = _series(case, "inletPressure")
    mom_i = _series(case, "inletMomentum")
    force = _series(case, "wallForce", "force.dat")

    n = min(len(a) for a in (flux_e, pres_e, mom_e, flux_i, pres_i, mom_i))
    it = flux_e[:n, 0]
    area_e, area_i = flux_e[:n, 1], flux_i[:n, 1]
    k_e = math.pi*re**2/area_e
    k_i = math.pi*rc**2/area_i

    mdot_e = flux_e[:n, 2]*k_e
    mdot_i = -flux_i[:n, 2]*k_i
    f_exit = (mom_e[:n, 1] + pres_e[:n, 1] - p_amb*area_e)*k_e
    # Inflow: phi < 0, so sum(phi U_x) is minus the momentum flux.
    f_inlet = (-mom_i[:n, 1] + pres_i[:n, 1] - p_amb*area_i)*k_i
    m = min(n, len(force))
    f_wall = np.full(n, np.nan)
    f_wall[:m] = force[:m, 1]*k_e[:m]   # total force, x component
    f_budget = f_inlet - f_wall

    last = it >= it[-1]*(1 - window)

    def stat(a):
        a = a[last]
        mean = float(np.mean(a))
        return mean, float((np.max(a) - np.min(a))/abs(mean))

    mdot, mdot_var = stat(mdot_e)
    thrust, thrust_var = stat(f_exit)
    mdot_in, _ = stat(mdot_i)
    budget, _ = stat(f_budget)
    j = min(len(avg_e), n) - 1
    design = rec["design"]
    try:
        thr = _series(case, "throatFlux")   # Time, area, areaNormalIntegrate
        rt = rec["geometry"]["rt"]
        mdot_t = thr[:, 2]*math.pi*rt**2/thr[:, 1]
        sel = thr[:, 0] >= thr[-1, 0]*(1 - window)
        mdot_throat = float(np.mean(mdot_t[sel]))
    except FileNotFoundError:
        mdot_throat = float("nan")
    result = {
        "iterations": float(it[-1]),
        "mdot_exit": mdot,
        "mdot_inlet": mdot_in,
        "mdot_throat": mdot_throat,
        "mdot_balance": (mdot_in - mdot)/mdot,
        "thrust_exit": thrust,
        "thrust_wall_budget": budget,
        "momentum_closure": (budget - thrust)/thrust,
        "isp": thrust/(mdot*G0),
        "pe_mean": float(avg_e[j, 1]),
        "te_mean": float(avg_e[j, 2]),
        "mach_e_mean": float(avg_e[j, 3]),
        "steadiness_mdot": mdot_var,
        "steadiness_thrust": thrust_var,
        # Injector-face mass flow: it settles first; equal to the exit and
        # throat values once mdot_balance is small.
        "discharge_coefficient": mdot_in/rec["mdot_ideal"],
        "design": {
            "model": design["model"],
            "thrust": design["thrust"],
            "mdot": design["mdot"],
            "isp": design["isp"],
            "pe": design["pe"],
            "te": design["te"],
            "mach_e": design["mach_e"],
        },
    }
    if not rec["explicit_geometry"]:
        result["thrust_error"] = (thrust - design["thrust"])/design["thrust"]
        result["mdot_error"] = (mdot - design["mdot"])/design["mdot"]
        result["pe_error"] = (result["pe_mean"] - design["pe"])/design["pe"]
    return result


def write(case, window=0.1):
    result = evaluate(case, window)
    with open(os.path.join(case, "results.json"), "w") as f:
        json.dump(result, f, indent=2)
    return result
