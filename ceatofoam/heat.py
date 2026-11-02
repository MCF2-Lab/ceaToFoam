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

"""Gas-side wall heat transfer: the Bartz correlation along a contour.

Equations of Huzel and Huang (1992), Modern Engineering for Design of
Liquid-Propellant Rocket Engines, pp. 85-86, in SI units:

* (4-10)   q = h_g (T_aw - T_wg)
* (4-10-a) T_aw = T_c [1 + r (gamma-1)/2 M^2] / [1 + (gamma-1)/2 M^2]
* r = Pr^0.5 laminar, r = Pr^0.33 turbulent
* (4-13)   h_g = [0.026/D_t^0.2 (mu^0.2 c_p/Pr^0.6) (p_c/c*)^0.8
                  (D_t/R)^0.1] (A_t/A)^0.9 sigma
* (4-14)   sigma, the correction for property variation across the layer
* (4-15)   Pr = 4 gamma/(9 gamma - 5), when no transport data are available

The properties in (4-13) are those of the chamber. Which gas supplies them
is the choice this module exposes: the frozen and the equilibrium c_p of a
rocket chamber differ by a factor between two and four, since the
equilibrium value carries the heat of the reactions.
"""

from dataclasses import asdict, dataclass

import numpy as np

from . import sizing


@dataclass
class ChamberGas:
    """Chamber properties for the correlation (SI)."""

    label: str
    mu: float
    cp: float
    prandtl: float
    gamma: float


def frozen_properties(mixture, tc):
    """Properties of the frozen mixture at the chamber temperature."""
    mixture._set(tc)
    gas = mixture.gas
    return ChamberGas(
        "frozen", gas.viscosity, gas.cp_mass,
        gas.viscosity*gas.cp_mass/gas.thermal_conductivity, mixture.gamma(tc),
    )


def equilibrium_properties(cea_point):
    """Properties CEA reports for the equilibrium chamber."""
    return ChamberGas(
        "equilibrium", cea_point.visc_c, cea_point.cp_c, cea_point.prandtl_c,
        cea_point.gammas_c,
    )


def prandtl_estimate(gamma):
    """Eq. (4-15)."""
    return 4*gamma/(9*gamma - 5)


def recovery_factor(prandtl, turbulent=True):
    """Local recovery factor (Huzel and Huang, p. 85)."""
    return prandtl**(0.33 if turbulent else 0.5)


def adiabatic_wall_temperature(tc, gamma, mach, recovery):
    """Eq. (4-10-a)."""
    m = 0.5*(gamma - 1)*mach**2
    return tc*(1 + recovery*m)/(1 + m)


def sigma(t_wall, tc, gamma, mach):
    """Eq. (4-14)."""
    m = 0.5*(gamma - 1)*mach**2
    return 1/((0.5*t_wall/tc*(1 + m) + 0.5)**0.68*(1 + m)**0.12)


def bartz(gas, pc, cstar, throat_diameter, throat_curvature, area_ratio,
          sigma_):
    """Gas-side coefficient h_g [W/(m2 K)], Eq. (4-13), at A/A_t."""
    return (
        0.026/throat_diameter**0.2
        *(gas.mu**0.2*gas.cp/gas.prandtl**0.6)
        *(pc/cstar)**0.8
        *(throat_diameter/throat_curvature)**0.1
    )*area_ratio**-0.9*sigma_


def along_contour(contour, flow_gas, tc, pc, cstar, gases, t_wall, x):
    """Bartz along the wall at axial positions ``x``.

    ``flow_gas`` gives the Mach number (1-D, frozen interface); ``gases`` is
    a list of ChamberGas, each giving h_g, T_aw and q; ``t_wall`` is a
    number or an array at ``x``. Returns a dict of arrays.
    """
    x = np.asarray(x, dtype=float)
    t_w = np.broadcast_to(np.asarray(t_wall, dtype=float), x.shape)
    r = contour.radius(x)
    throat = sizing.frozen_throat(flow_gas, tc, pc)
    mach = np.empty_like(x)
    for k, (xi, ri) in enumerate(zip(x, r)):
        mach[k] = sizing.frozen_station_at_area(
            flow_gas, tc, pc, float((ri/contour.rt)**2),
            xi > contour.x_throat, throat,
        ).mach
    out = {"x": x, "r": r, "mach": mach, "t_wall": np.array(t_w)}
    # Throat curvature: the mean of the arcs on either side.
    curvature = 0.5*(contour.r_tu + contour.r_td)
    for g in gases:
        rec = recovery_factor(g.prandtl)
        t_aw = adiabatic_wall_temperature(tc, g.gamma, mach, rec)
        sig = sigma(t_w, tc, g.gamma, mach)
        h = bartz(g, pc, cstar, 2*contour.rt, curvature,
                  (r/contour.rt)**2, sig)
        out[f"t_aw_{g.label}"] = t_aw
        out[f"h_{g.label}"] = h
        out[f"q_{g.label}"] = h*(t_aw - t_w)
    return out


def record(gases):
    """Input record of the correlation."""
    return {g.label: asdict(g) for g in gases}
