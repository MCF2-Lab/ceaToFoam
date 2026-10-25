# ceaToFoam: NASA CEA chamber thermochemistry to OpenFOAM rocket-nozzle cases
# Copyright (C) 2025-2026 Luke Cerda and Ehsan Madadi-Kandjani
# SPDX-License-Identifier: GPL-3.0-only
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License version 3 as published
# by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

"""Thrust, mass flow and exit state from a finished case.

Wedge integrals are scaled to the full annulus with the ratio of the
circular area pi r^2 to the measured patch area (the same ratio holds at
every radius of a planar-faced wedge).

Thrust is evaluated in two ways:

* exit plane: F = sum(phi U_x) + integral (p - p_amb) dA over the exit faceZone;
* wall budget: F = [sum(phi U_x) + integral (p - p_amb) dA] at the injector
  face minus the axial force of the gas on the wall.

Both are equal if axial momentum is conserved, so their difference measures
the momentum budget closure of the discretization. Results are means over
the last ``window`` fraction of the iterations; their relative range
(steadiness) and relative standard deviation (iterative uncertainty) are
reported with them.
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


def _read_minmax_T(case):
    """(iteration, min T, max T) rows from the minMax function object."""
    base = os.path.join(case, "postProcessing", "minMax")
    rows = []
    for t in sorted(os.listdir(base), key=float):
        path = os.path.join(base, t, "fieldMinMax.dat")
        if not os.path.exists(path):
            continue
        with open(path) as f:
            for line in f:
                tok = line.split()
                if len(tok) > 6 and not line.startswith("#") and tok[1] == "T":
                    # Time field min (location) max (location)
                    vals = line.replace("(", " ").replace(")", " ").split()
                    rows.append((float(vals[0]), float(vals[2]),
                                 float(vals[6])))
    return np.array(rows)


def evaluate(case, window=0.25):
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
        """Mean, relative range and relative standard deviation."""
        a = a[last]
        mean = float(np.mean(a))
        return (mean, float((np.max(a) - np.min(a))/abs(mean)),
                float(np.std(a)/abs(mean)))

    mdot, mdot_var, mdot_std = stat(mdot_e)
    thrust, thrust_var, thrust_std = stat(f_exit)
    mdot_in, _, _ = stat(mdot_i)
    budget, _, _ = stat(f_budget)
    j = min(len(avg_e), n) - 1
    design = rec["design"]
    # Extremes of T in the averaging window: start-up transients may touch
    # the janaf limits (clipped by OpenFOAM); the converged state must not.
    try:
        mm = _read_minmax_T(case)
        sel = mm[:, 0] >= mm[-1, 0]*(1 - window)
        t_min_window = float(mm[sel, 1].min())
        t_max_window = float(mm[sel, 2].max())
    except (FileNotFoundError, IndexError, ValueError):
        t_min_window = t_max_window = float("nan")
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
        # Iterative uncertainty: relative standard deviation in the window.
        "iterative_std_mdot": mdot_std,
        "iterative_std_thrust": thrust_std,
        # Injector-face mass flow: it settles first; equal to the exit and
        # throat values once mdot_balance is small.
        "discharge_coefficient": mdot_in/rec["mdot_ideal"],
        "t_min_window": t_min_window,
        "t_max_window": t_max_window,
        "janaf_t_low": rec["fits"]["janaf"].get("t_low", 200.0),
        "design": {
            "model": design["model"],
            "thrust": design["thrust"],
            "mdot": design["mdot"],
            "isp": design["isp"],
            "pe": design["pe"],
            "te": design["te"],
            "mach_e": design["mach_e"],
            "gamma_c": design["gamma_c"],
            "area_ratio": design["area_ratio"],
            "cstar": design["cstar"],
        },
    }
    if not rec["explicit_geometry"]:
        result["thrust_error"] = (thrust - design["thrust"])/design["thrust"]
        result["mdot_error"] = (mdot - design["mdot"])/design["mdot"]
        result["pe_error"] = (result["pe_mean"] - design["pe"])/design["pe"]
    return result


def profiles(case):
    """Latest wall and axis profiles as CSV text: x, x - x_throat, r,
    p/p0, T, Ma (p0 the chamber pressure). Returns {name: csv_text}."""
    with open(os.path.join(case, "ceatofoam.json")) as f:
        rec = json.load(f)
    base = os.path.join(case, "postProcessing", "profiles")
    if not os.path.isdir(base):
        return {}
    latest = max(os.listdir(base), key=float)
    p0 = rec["chamber"]["p"]
    x_t = rec["geometry"]["x_throat"]
    out = {}
    for name in ("wall", "axis"):
        path = os.path.join(base, latest, f"{name}_Ma_T_p_U.xy")
        if not os.path.exists(path):
            continue
        data = np.loadtxt(path)
        x = data[:, 0]
        if name == "wall":
            stations = np.array(rec["geometry"]["stations"])
            r = np.interp(x, stations[:, 0], stations[:, 1])
        else:
            r = np.zeros_like(x)
        rows = ["x,x_minus_throat,r,p_over_p0,T,Ma"] + [
            f"{a:.8g},{a - x_t:.8g},{b:.8g},{c/p0:.8g},{d:.8g},{e:.8g}"
            for a, b, c, d, e in zip(x, r, data[:, 3], data[:, 2],
                                     data[:, 1])
        ]
        out[name] = "\n".join(rows) + "\n"
    return out


def write(case, window=0.25):
    result = evaluate(case, window)
    for name, text in profiles(case).items():
        with open(os.path.join(case, f"profile_{name}.csv"), "w") as f:
            f.write(text)
    with open(os.path.join(case, "results.json"), "w") as f:
        json.dump(result, f, indent=2)
    return result
