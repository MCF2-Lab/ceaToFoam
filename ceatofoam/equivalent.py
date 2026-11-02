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

"""Equivalent gas: a thermally perfect gas with the expansion of a reacting one.

Along a streamline of an inviscid adiabatic flow the first law reads
de = (p/rho^2) drho, with e the total internal energy per unit mass, heat of
formation included, whatever the chemistry does. Define for a fixed gas
constant R the temperature tau = p/(rho R) and the function e_eq(tau) equal
to e along a given expansion path p(rho), e(rho). The thermally perfect gas
p = rho R tau, e = e_eq(tau) then has the same first law along the same
path, so its isentrope through the chamber is that path, and its speed of
sound is dp/drho along it. It is a gas that rhoCentralFoam, janaf and
perfectGas solve with no change.

Two paths are available: shifting equilibrium (Cantera at constant entropy)
and a kinetic expansion through a given contour (Cantera rates; the marching
variable is the pressure, so the sonic point is not singular). The gas
constant is that of the chamber, so tau = T at the chamber.

tau is not the gas temperature, which the path gives as a function of tau
(EquivalentGas.true_temperature). The gas is exact on the path and only
there: it does not hold across strong shocks or inside wall layers, where
the state leaves the path, and it has one path per chamber state, as the
frozen mixture has one composition.

EquivalentGas has the interface of thermo.FrozenMixture, so the frozen
sizing functions evaluate its 1-D expansion and reproduce the path; this is
the check of the fit that build.make_case records.
"""

from dataclasses import dataclass, field
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .thermo import JanafFit
from .units import P_STD, RU

#: Fractions of the tau range tried for the common temperature of the fit.
_COMMON_FRACTIONS = (0.3, 0.4, 0.5, 0.6)
#: Weight of the constant continuation of cp outside the path, which keeps
#: the quartics finite in the start-up range without loading the fit.
_TAIL_WEIGHT = 0.05
#: Range of the fit relative to the path: start-up transients overshoot the
#: steady extremes of tau (test case: 10 % above T_c).
_TAU_MIN_FACTOR = 0.7
_TAU_MAX_FACTOR = 1.15


@dataclass
class ExpansionPath:
    """States of one expansion from the chamber, in order of falling p."""

    kind: str                   # equilibrium | kinetic
    p: np.ndarray               # [Pa]
    rho: np.ndarray             # [kg/m3]
    e: np.ndarray               # total internal energy [J/kg]
    h: np.ndarray               # total enthalpy [J/kg]
    T: np.ndarray               # gas temperature [K]
    mw: np.ndarray              # molecular weight [kg/kmol]
    mu: np.ndarray              # viscosity [Pa s]
    extra: dict = field(default_factory=dict)

    @property
    def velocity(self):
        return np.sqrt(np.maximum(2*(self.h[0] - self.h), 0.0))

    @property
    def mass_flux(self):
        return self.rho*self.velocity


def _gas(composition, tc, pc, mechanism):
    """Cantera gas in equilibrium at (tc, pc)."""
    import cantera as ct

    gas = ct.Solution(mechanism)
    gas.TPX = tc, pc, composition
    gas.equilibrate("TP")
    return gas


def equilibrium_path(composition, tc, pc, p_end, mechanism="gri30_highT.yaml",
                     n=1200):
    """Shifting-equilibrium expansion at the chamber entropy, pc to p_end."""
    gas = _gas(composition, tc, pc, mechanism)
    s0 = gas.entropy_mass
    rows = []
    for p in np.geomspace(pc, p_end, n):
        gas.SP = s0, p
        gas.equilibrate("SP")
        rows.append((p, gas.density, gas.int_energy_mass, gas.enthalpy_mass,
                     gas.T, gas.mean_molecular_weight, gas.viscosity))
    a = np.array(rows).T
    return ExpansionPath("equilibrium", *a)


def _axial_position(path, contour):
    """Axial position of each state of a path in ``contour`` [m]: the
    area ratio from the mass flux, the wall radius on either side of the
    throat; beyond the exit the divergent cone is continued."""
    g = path.mass_flux
    g[0] = max(g[0], 1e-9)
    i_t = int(np.argmax(g))
    r = contour.rt*np.sqrt(g[i_t]/np.maximum(g, 1e-9))
    x_grid = np.linspace(0.0, contour.x_exit, 2001)
    r_grid = contour.radius(x_grid)
    j_t = int(np.argmin(np.abs(x_grid - contour.x_throat)))
    x = np.empty_like(r)
    x[:i_t + 1] = np.interp(-r[:i_t + 1], -r_grid[:j_t + 1], x_grid[:j_t + 1])
    x[i_t + 1:] = np.interp(r[i_t + 1:], r_grid[j_t:], x_grid[j_t:])
    beyond = r > contour.re
    x[beyond & (np.arange(len(r)) > i_t)] = contour.x_exit + (
        r[beyond & (np.arange(len(r)) > i_t)] - contour.re
    )/math.tan(math.radians(contour.theta_d))
    return x


def _march(gas, h0, xi, x_of_xi, T0, Y0):
    """Integrate T and Y against xi = -ln p for a prescribed x(xi)."""
    dx_dxi = np.gradient(x_of_xi, xi)

    def rhs(l, y):
        p = math.exp(-l)
        gas.TPY = y[0], p, y[1:]
        u = math.sqrt(max(2*(h0 - gas.enthalpy_mass), 1.0))
        dt_dxi = max(float(np.interp(l, xi, dx_dxi)), 0.0)/u
        rate = gas.net_production_rates*gas.molecular_weights/gas.density
        hi = gas.partial_molar_enthalpies/gas.molecular_weights
        # cp dT + sum h_i dY_i = dp/rho, with dp = -p dxi
        dT = (-p/gas.density - float(np.dot(hi, rate))*dt_dxi)/gas.cp_mass
        return np.concatenate(([dT], rate*dt_dxi))

    sol = solve_ivp(
        rhs, (xi[0], xi[-1]), np.concatenate(([T0], Y0)), method="BDF",
        t_eval=xi, rtol=1e-7, atol=1e-11,
    )
    if not sol.success:
        raise RuntimeError(f"kinetic march failed: {sol.message}")
    return sol.y


def kinetic_path(composition, tc, pc, p_end, contour,
                 mechanism="gri30_highT.yaml", n=400, passes=2):
    """Finite-rate expansion through ``contour`` from the inlet station.

    Shifting equilibrium holds in the chamber, up to the pressure at which
    the path reaches the contraction ratio of the contour. From there the
    composition is integrated with the rates of ``mechanism``, with the
    pressure as marching variable and x(p) from the area of the previous
    pass (first pass: the equilibrium path), so that the residence time
    belongs to the contour and the sonic point is an ordinary station.
    """
    gas = _gas(composition, tc, pc, mechanism)
    if not gas.n_reactions:
        raise ValueError(f"{mechanism} has no reactions")
    h0 = gas.enthalpy_mass
    eq = equilibrium_path(composition, tc, pc, p_end, mechanism)
    g = eq.mass_flux
    i_t = int(np.argmax(g))
    ratio = g[i_t]/np.maximum(g[:i_t + 1], 1e-9)
    i0 = int(np.argmin(np.abs(ratio - contour.contraction_ratio)))

    # Start state: equilibrium at the inlet pressure on the chamber isentrope.
    s0 = gas.entropy_mass
    gas.SP = s0, eq.p[i0]
    gas.equilibrate("SP")
    T_start, Y_start = gas.T, gas.Y.copy()

    xi = np.linspace(-math.log(eq.p[i0]), -math.log(p_end), n)
    x_eq = _axial_position(eq, contour)
    # x against xi along the equilibrium path (interp needs rising xi)
    x_of_xi = np.interp(xi, -np.log(eq.p), x_eq)
    x_start = float(x_of_xi[0])
    x_of_xi = np.maximum.accumulate(x_of_xi)

    prev = None
    for _ in range(max(1, passes)):
        y = _march(gas, h0, xi, x_of_xi, T_start, Y_start)
        rows = []
        for k, l in enumerate(xi):
            p = math.exp(-l)
            gas.TPY = y[0, k], p, y[1:, k]
            rows.append((p, gas.density, gas.int_energy_mass,
                         gas.enthalpy_mass, gas.T,
                         gas.mean_molecular_weight, gas.viscosity))
        a = np.array(rows).T
        head = eq.p > eq.p[i0]*(1 + 1e-9)
        merged = [np.concatenate((b[head], c)) for b, c in zip(
            (eq.p, eq.rho, eq.e, eq.h, eq.T, eq.mw, eq.mu), a)]
        path = ExpansionPath("kinetic", *merged)
        x_new = _axial_position(path, contour)
        x_of_xi = np.maximum.accumulate(np.maximum(
            np.interp(xi, -np.log(path.p), x_new), x_start
        ))
        prev = x_of_xi
    path.extra = {"x": np.interp(-np.log(path.p), xi, prev).tolist(),
                  "residence_time": _residence_time(path, xi, prev)}
    return path


def _residence_time(path, xi, x_of_xi):
    """Time from the inlet station to the end of the path [s]."""
    u = np.interp(xi, -np.log(path.p), path.velocity)
    return float(np.trapezoid(1/np.maximum(u, 1.0), x_of_xi))


def fit_path(tau, cp, R, t_low, t_high):
    """Two janaf quartics to cp(tau), continuous with their slope, on the
    path; below and above it a low-weight constant continuation."""
    y = cp/R
    tails_lo = np.arange(t_low, tau[0], 10.0)
    tails_hi = np.arange(tau[-1] + 10.0, t_high + 5.0, 10.0)
    T = np.concatenate((tails_lo, tau, tails_hi))
    Y = np.concatenate((np.full(len(tails_lo), y[0]), y,
                        np.full(len(tails_hi), y[-1])))
    W = np.concatenate((np.full(len(tails_lo), _TAIL_WEIGHT),
                        np.ones(len(tau)),
                        np.full(len(tails_hi), _TAIL_WEIGHT)))
    basis = lambda t: np.vander(np.atleast_1d(t), 5, increasing=True)
    dbasis = lambda t: np.array([[0.0, 1.0, 2*t, 3*t**2, 4*t**3]])

    def solve(t_common):
        lo = T < t_common
        a = np.zeros((len(T), 10))
        a[lo, :5] = basis(T[lo])
        a[~lo, 5:] = basis(T[~lo])
        c = np.zeros((2, 10))
        c[0, :5], c[0, 5:] = basis(t_common)[0], -basis(t_common)[0]
        c[1, :5], c[1, 5:] = dbasis(t_common)[0], -dbasis(t_common)[0]
        scale = np.concatenate([t_high**-np.arange(5.0)]*2)
        a_s, c_s = (a*scale)*W[:, None], c*scale
        kkt = np.block([[2*a_s.T @ a_s, c_s.T], [c_s, np.zeros((2, 2))]])
        rhs = np.concatenate((2*a_s.T @ (Y*W), np.zeros(2)))
        x = np.linalg.solve(kkt, rhs)[:10]*scale
        low, high = np.zeros(7), np.zeros(7)
        low[:5], high[:5] = x[:5], x[5:]
        high[5] = (JanafFit.h_over_r(low, t_common)
                   - JanafFit.h_over_r(high, t_common))
        high[6] = (JanafFit.s_over_r(low, t_common)
                   - JanafFit.s_over_r(high, t_common))
        fit_y = np.array([
            JanafFit.cp_over_r(low if t < t_common else high, t) for t in tau
        ])
        return float(np.max(np.abs(fit_y/y - 1))), low, high

    best = None
    for f in _COMMON_FRACTIONS:
        t_common = tau[0] + f*(tau[-1] - tau[0])
        err, low, high = solve(t_common)
        if best is None or err < best[0]:
            best = (err, low, high, t_common)
    err, low, high, t_common = best
    return JanafFit(
        t_low=t_low, t_common=t_common, t_high=t_high, low=low, high=high,
        cp_max_rel_error=err, h_max_abs_error=float("nan"),
    )


class EquivalentGas:
    """Thermally perfect gas with the expansion ``path`` as its isentrope.

    Interface of thermo.FrozenMixture, with temperature standing for tau.
    """

    def __init__(self, path, mixture, tc):
        self.path = path
        self.kind = path.kind
        self.mw = float(path.mw[0])
        self.R = RU/self.mw
        tau = path.p/(path.rho*self.R)
        if np.any(np.diff(tau) >= 0.0):
            raise ValueError("tau does not fall along the path")
        order = np.argsort(tau)
        self.tau = tau[order]
        e = path.e[order]
        self.cv_path = np.gradient(e, self.tau)
        if np.any(self.cv_path <= 0.0):
            raise ValueError("cv of the equivalent gas is not positive")
        self.cp_path = self.cv_path + self.R
        self.true_t = path.T[order]
        self.true_mw = path.mw[order]
        self.mu_path = path.mu[order]
        # Trim the end points, where the one-sided differences sit.
        keep = slice(2, -2)
        self.janaf = fit_path(
            self.tau[keep], self.cp_path[keep], self.R,
            _TAU_MIN_FACTOR*self.tau[0], _TAU_MAX_FACTOR*tc,
        )
        self.gas_tc = tc
        self.mechanism = mixture.mechanism
        self.composition = mixture.composition
        self.dropped = mixture.dropped
        self.dropped_fraction = mixture.dropped_fraction
        self.t_data_min = float(self.janaf.t_low)

    def _a(self, T):
        return self.janaf.low if T < self.janaf.t_common else self.janaf.high

    def cp(self, T):
        return self.R*JanafFit.cp_over_r(self._a(T), T)

    def h(self, T):
        return self.R*JanafFit.h_over_r(self._a(T), T)

    def s(self, T, p):
        return self.R*(JanafFit.s_over_r(self._a(T), T) - math.log(p/P_STD))

    def gamma(self, T):
        cp = self.cp(T)
        return cp/(cp - self.R)

    def mu(self, T):
        return float(np.interp(T, self.tau, self.mu_path))

    def true_temperature(self, tau):
        """Gas temperature [K] at ``tau`` on the path."""
        return np.interp(tau, self.tau, self.true_t)

    def temperature_at(self, s0, p, bracket=None):
        lo, hi = bracket or (self.janaf.t_low, self.janaf.t_high)
        return brentq(lambda T: self.s(T, p) - s0, lo, hi, xtol=1e-10)

    def record(self, n=41):
        """Summary for ceatofoam.json, with the map from tau to T."""
        grid = np.linspace(self.tau[0], self.tau[-1], n)
        path = self.path
        return {
            "path": self.kind,
            "molecular_weight": self.mw,
            "tau": grid.tolist(),
            "temperature": np.interp(grid, self.tau, self.true_t).tolist(),
            "path_molecular_weight":
                np.interp(grid, self.tau, self.true_mw).tolist(),
            "gamma": [self.gamma(float(t)) for t in grid],
            "cp_max_rel_error": self.janaf.cp_max_rel_error,
            "residence_time": path.extra.get("residence_time"),
        }


def build(kind, mixture, tc, pc, p_end, contour=None):
    """EquivalentGas of the chamber state: ``kind`` equilibrium or kinetic.

    The path ends well below ``p_end``, so that nozzles larger than the one
    sized for it stay inside the range of the fit.
    """
    if kind == "equilibrium":
        path = equilibrium_path(
            mixture.composition, tc, pc, 0.15*p_end, mixture.mechanism
        )
    elif kind == "kinetic":
        if contour is None:
            raise ValueError("a kinetic path needs the contour")
        path = kinetic_path(
            mixture.composition, tc, pc, 0.15*p_end, contour,
            mixture.mechanism,
        )
    else:
        raise ValueError(f"unknown gas model {kind!r}")
    return EquivalentGas(path, mixture, tc)
