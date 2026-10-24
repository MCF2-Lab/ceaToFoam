"""Frozen chamber mixture: Cantera properties and OpenFOAM coefficient fits.

The CEA chamber mole fractions are frozen and evaluated with a Cantera
mechanism (default gri30_highT.yaml, NASA 7-coefficient thermodynamics
valid 200-6000 K, mixture-averaged transport). Species that the mechanism
does not contain are dropped and the remaining mole fractions renormalized;
the dropped mole fraction is reported.

Every quantity used for sizing and every coefficient written for OpenFOAM
comes from this one mixture, so the gas constant, cp(T) and gamma(T) of the
1-D design and of the CFD are the same.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .units import P_STD, RU, T_REF


class FrozenMixture:
    """Ideal-gas mixture of fixed composition."""

    def __init__(self, composition, mechanism="gri30_highT.yaml"):
        import cantera as ct

        self.gas = ct.Solution(mechanism)
        names = set(self.gas.species_names)
        kept = {s: x for s, x in composition.items() if s in names and x > 0}
        if not kept:
            raise ValueError(f"no species of {composition} in {mechanism}")
        total = sum(kept.values())
        self.dropped = {
            s: x for s, x in composition.items() if s not in names and x > 0
        }
        self.dropped_fraction = sum(self.dropped.values())/(
            total + sum(self.dropped.values())
        )
        self.composition = {s: x/total for s, x in kept.items()}
        self.mechanism = mechanism
        #: Lowest temperature covered by the thermodynamic data of every
        #: species in the mixture; below it the polynomials are extrapolated.
        self.t_data_min = max(
            self.gas.species(s).thermo.min_temp for s in self.composition
        )
        self.gas.TPX = 300.0, P_STD, self.composition
        #: Molecular weight [kg/kmol] of the frozen mixture.
        self.mw = self.gas.mean_molecular_weight
        #: Specific gas constant [J/(kg K)].
        self.R = RU/self.mw

    def _set(self, T, p=P_STD):
        self.gas.TP = float(T), float(p)

    def cp(self, T):
        """Frozen specific heat at constant pressure [J/(kg K)]."""
        self._set(T)
        return self.gas.cp_mass

    def h(self, T):
        """Absolute specific enthalpy [J/kg] (ideal gas: independent of p)."""
        self._set(T)
        return self.gas.enthalpy_mass

    def s(self, T, p):
        """Specific entropy [J/(kg K)] at temperature T and pressure p."""
        self._set(T, p)
        return self.gas.entropy_mass

    def gamma(self, T):
        """Frozen ratio of specific heats cp/(cp - R)."""
        cp = self.cp(T)
        return cp/(cp - self.R)

    def mu(self, T):
        """Mixture-averaged dynamic viscosity [Pa s]."""
        self._set(T)
        return self.gas.viscosity

    def temperature_at(self, s0, p, bracket=(50.0, 6000.0)):
        """Temperature where s(T, p) = s0 (isentropic state at pressure p)."""
        return brentq(lambda T: self.s(T, p) - s0, *bracket, xtol=1e-10)


@dataclass
class JanafFit:
    """Two-range JANAF coefficients in OpenFOAM's convention.

    ``low`` and ``high`` hold (a0..a4, a5, a6) with Cp/R a quartic in T,
    H/R = sum a_k T^(k+1)/(k+1) + a5 and S/R = a0 ln T + ... + a6 at the
    standard pressure; OpenFOAM multiplies them by its R = RU/molWeight.
    """

    t_low: float
    t_common: float
    t_high: float
    low: np.ndarray
    high: np.ndarray
    cp_max_rel_error: float
    h_max_abs_error: float

    @staticmethod
    def cp_over_r(a, T):
        return (((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]

    @staticmethod
    def h_over_r(a, T):
        return (
            (((a[4]/5*T + a[3]/4)*T + a[2]/3)*T + a[1]/2)*T + a[0]
        )*T + a[5]

    @staticmethod
    def s_over_r(a, T):
        return (
            (((a[4]/4*T + a[3]/3)*T + a[2]/2)*T + a[1])*T + a[0]*np.log(T)
            + a[6]
        )

    def coefficients(self, T):
        return self.low if T < self.t_common else self.high


def fit_janaf(mixture, t_low=200.0, t_common=1000.0, t_high=6000.0, step=5.0):
    """Fit OpenFOAM janaf coefficients to the frozen mixture.

    The two Cp/R quartics are fitted by least squares with the constraint
    that Cp/R and its slope are continuous at ``t_common``. The low-range
    enthalpy and entropy constants reproduce the mixture's h and s at
    T_REF and the standard pressure; the high-range constants make h and s
    continuous at ``t_common``. The fit errors over [t_low, t_high] are
    returned with the coefficients.
    """
    R = mixture.R
    t1 = np.arange(t_low, t_common + 0.5*step, step)
    t2 = np.arange(t_common, t_high + 0.5*step, step)
    y1 = np.array([mixture.cp(T) for T in t1])/R
    y2 = np.array([mixture.cp(T) for T in t2])/R

    def basis(T):
        return np.vander(np.atleast_1d(T), 5, increasing=True)

    def dbasis(T):
        return np.array([[0.0, 1.0, 2*T, 3*T**2, 4*T**3]])

    # Unknowns: 5 low + 5 high coefficients. Equality constraints C x = 0.
    a_mat = np.zeros((len(t1) + len(t2), 10))
    a_mat[:len(t1), :5] = basis(t1)
    a_mat[len(t1):, 5:] = basis(t2)
    rhs = np.concatenate([y1, y2])
    c_mat = np.zeros((2, 10))
    c_mat[0, :5], c_mat[0, 5:] = basis(t_common)[0], -basis(t_common)[0]
    c_mat[1, :5], c_mat[1, 5:] = dbasis(t_common)[0], -dbasis(t_common)[0]
    # Scale columns for conditioning, solve the KKT system.
    scale = np.concatenate([t_high**-np.arange(5.0)]*2)
    a_s, c_s = a_mat*scale, c_mat*scale
    kkt = np.block([
        [2*a_s.T @ a_s, c_s.T],
        [c_s, np.zeros((2, 2))],
    ])
    sol = np.linalg.solve(
        kkt, np.concatenate([2*a_s.T @ rhs, np.zeros(2)])
    )
    x = sol[:10]*scale
    low = np.zeros(7)
    high = np.zeros(7)
    low[:5], high[:5] = x[:5], x[5:]

    low[5] = mixture.h(T_REF)/R - JanafFit.h_over_r(low, T_REF)
    low[6] = mixture.s(T_REF, P_STD)/R - JanafFit.s_over_r(low, T_REF)
    high[5] = (
        JanafFit.h_over_r(low, t_common) - JanafFit.h_over_r(high, t_common)
    )
    high[6] = (
        JanafFit.s_over_r(low, t_common) - JanafFit.s_over_r(high, t_common)
    )

    t_all = np.concatenate([t1, t2[1:]])
    cp_true = np.concatenate([y1, y2[1:]])
    cp_fit = np.array([
        JanafFit.cp_over_r(low if T < t_common else high, T) for T in t_all
    ])
    h_true = np.array([mixture.h(T) for T in t_all])/R
    h_fit = np.array([
        JanafFit.h_over_r(low if T < t_common else high, T) for T in t_all
    ])
    return JanafFit(
        t_low=t_low,
        t_common=t_common,
        t_high=t_high,
        low=low,
        high=high,
        cp_max_rel_error=float(np.max(np.abs(cp_fit/cp_true - 1))),
        h_max_abs_error=float(np.max(np.abs(h_fit - h_true))*R),
    )


@dataclass
class SutherlandFit:
    """mu = As sqrt(T)/(1 + Ts/T), fitted to the mixture viscosity."""

    As: float
    Ts: float
    t_min: float
    t_max: float
    max_rel_error: float

    def mu(self, T):
        return self.As*np.sqrt(T)/(1 + self.Ts/T)


def fit_sutherland(mixture, t_min, t_max, n=200):
    """Fit Sutherland coefficients over [t_min, t_max].

    sqrt(T)/mu = 1/As + (Ts/As)/T is linear in 1/T; it is fitted by least
    squares weighted by mu so that the residual is the relative error of mu.
    """
    T = np.linspace(t_min, t_max, n)
    mu = np.array([mixture.mu(t) for t in T])
    y = np.sqrt(T)/mu
    a_mat = np.column_stack([np.ones_like(T), 1/T])
    w = 1/y
    b, *_ = np.linalg.lstsq(a_mat*w[:, None], y*w, rcond=None)
    As, Ts = 1/b[0], b[1]/b[0]
    fit = As*np.sqrt(T)/(1 + Ts/T)
    return SutherlandFit(
        As=float(As),
        Ts=float(Ts),
        t_min=float(t_min),
        t_max=float(t_max),
        max_rel_error=float(np.max(np.abs(fit/mu - 1))),
    )
