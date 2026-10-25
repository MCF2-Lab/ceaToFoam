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

"""Design-point sizing of the throat and exit.

Three sizing models are available. All give a nozzle that is matched
(p_e = p_amb) at the design point and deliver the design thrust there:

``frozen`` (default)
    Isentropic expansion of the frozen chamber mixture with temperature-
    dependent cp (thermo.FrozenMixture). The throat is the sonic station.
    This is the thermodynamics the CFD solves, so it is the consistent model.
``equilibrium``
    CEA shifting-equilibrium performance at pc/p_amb: exit velocity from Isp,
    throat from c*, area ratio from CEA.
``constant-gamma``
    The textbook constant-gamma relations (Sutton and Biblarz) with a
    chosen gamma and gas constant; gamma = CEA chamber gamma_s reproduces
    the sizing of ceaToFoam 1.x.
"""

from dataclasses import asdict, dataclass
import math

from scipy.optimize import brentq

from .units import G0


@dataclass
class Design:
    """Design-point quantities (SI)."""

    model: str
    thrust: float
    pc: float
    tc: float
    p_amb: float
    pe: float
    te: float
    mach_e: float
    ve: float
    mdot: float
    cstar: float
    area_ratio: float
    at: float
    ae: float
    gamma_c: float
    gamma_e: float
    mw: float

    @property
    def rt(self):
        return math.sqrt(self.at/math.pi)

    @property
    def re(self):
        return math.sqrt(self.ae/math.pi)

    @property
    def isp(self):
        return self.ve/G0

    def as_dict(self):
        d = asdict(self)
        d.update(rt=self.rt, re=self.re, isp=self.isp)
        return d


# ---------------------------------------------------------------------------
# Constant-gamma relations
# ---------------------------------------------------------------------------

def area_ratio_from_mach(mach, gamma):
    """A/A* of a calorically perfect gas at Mach number ``mach``."""
    g = gamma
    return (1/mach)*(
        (2/(g + 1))*(1 + 0.5*(g - 1)*mach**2)
    )**((g + 1)/(2*(g - 1)))


def mach_from_pressure_ratio(p0_over_p, gamma):
    """Mach number of an isentropic expansion to p = p0/p0_over_p."""
    g = gamma
    return math.sqrt((2/(g - 1))*(p0_over_p**((g - 1)/g) - 1))


def mass_flux_function(gamma):
    """Gamma(gamma) = sqrt(gamma) (2/(gamma+1))^((gamma+1)/(2(gamma-1)))."""
    g = gamma
    return math.sqrt(g)*(2/(g + 1))**((g + 1)/(2*(g - 1)))


def size_constant_gamma(thrust, pc, tc, p_amb, gamma, mw):
    from .units import RU

    R = RU/mw
    mach_e = mach_from_pressure_ratio(pc/p_amb, gamma)
    te = tc/(1 + 0.5*(gamma - 1)*mach_e**2)
    ve = math.sqrt(
        2*gamma*R*tc/(gamma - 1)*(1 - (p_amb/pc)**((gamma - 1)/gamma))
    )
    mdot = thrust/ve
    cstar = math.sqrt(R*tc)/mass_flux_function(gamma)
    at = mdot*cstar/pc
    eps = area_ratio_from_mach(mach_e, gamma)
    return Design(
        model="constant-gamma", thrust=thrust, pc=pc, tc=tc, p_amb=p_amb,
        pe=p_amb, te=te, mach_e=mach_e, ve=ve, mdot=mdot, cstar=cstar,
        area_ratio=eps, at=at, ae=eps*at, gamma_c=gamma, gamma_e=gamma,
        mw=mw,
    )


# ---------------------------------------------------------------------------
# Frozen, temperature-dependent cp
# ---------------------------------------------------------------------------

@dataclass
class Station:
    p: float
    T: float
    rho: float
    u: float
    a: float

    @property
    def mach(self):
        return self.u/self.a

    @property
    def mass_flux(self):
        return self.rho*self.u


def frozen_station(mixture, tc, pc, p):
    """Isentropic frozen state at pressure ``p`` from stagnation (tc, pc)."""
    s0 = mixture.s(tc, pc)
    h0 = mixture.h(tc)
    T = mixture.temperature_at(s0, p)
    u = math.sqrt(max(2*(h0 - mixture.h(T)), 0.0))
    return Station(
        p=p,
        T=T,
        rho=p/(mixture.R*T),
        u=u,
        a=math.sqrt(mixture.gamma(T)*mixture.R*T),
    )


def frozen_throat(mixture, tc, pc):
    """Sonic station of the frozen isentropic expansion."""
    p = brentq(
        lambda p: frozen_station(mixture, tc, pc, p).mach - 1,
        0.3*pc, 0.95*pc, xtol=1e-9*pc,
    )
    return frozen_station(mixture, tc, pc, p)


def frozen_station_at_area(mixture, tc, pc, area_ratio, supersonic,
                           throat=None):
    """Frozen isentropic state where A/A* = ``area_ratio``.

    ``supersonic`` selects the branch downstream of the throat.
    """
    throat = throat or frozen_throat(mixture, tc, pc)
    if area_ratio <= 1.0 + 1e-12:
        return throat
    target = throat.mass_flux/area_ratio
    f = lambda p: frozen_station(mixture, tc, pc, p).mass_flux - target
    if supersonic:
        # Mass flux falls monotonically with p on this branch.
        lo = 0.5*throat.p
        while f(lo) > 0.0:
            lo *= 0.5
        p = brentq(f, lo, throat.p, xtol=1e-12*pc)
    else:
        p = brentq(f, throat.p, pc*(1 - 1e-12), xtol=1e-12*pc)
    return frozen_station(mixture, tc, pc, p)


def size_frozen(thrust, pc, tc, p_amb, mixture):
    throat = frozen_throat(mixture, tc, pc)
    exit_ = frozen_station(mixture, tc, pc, p_amb)
    mdot = thrust/exit_.u
    at = mdot/throat.mass_flux
    eps = throat.mass_flux/exit_.mass_flux
    return Design(
        model="frozen", thrust=thrust, pc=pc, tc=tc, p_amb=p_amb,
        pe=p_amb, te=exit_.T, mach_e=exit_.mach, ve=exit_.u, mdot=mdot,
        cstar=pc/throat.mass_flux, area_ratio=eps, at=at, ae=eps*at,
        gamma_c=mixture.gamma(tc), gamma_e=mixture.gamma(exit_.T),
        mw=mixture.mw,
    )


def size_equilibrium(thrust, p_amb, cea_point):
    """Sizing from a CEA shifting-equilibrium point at pip = pc/p_amb."""
    pt = cea_point
    ve = pt.isp*G0
    mdot = thrust/ve
    at = mdot*pt.cstar/pt.pc
    return Design(
        model="equilibrium", thrust=thrust, pc=pt.pc, tc=pt.tc,
        p_amb=p_amb, pe=pt.p_e, te=pt.t_e, mach_e=pt.mach_e, ve=ve,
        mdot=mdot, cstar=pt.cstar, area_ratio=pt.area_ratio, at=at,
        ae=pt.area_ratio*at, gamma_c=pt.gammas_c, gamma_e=float("nan"),
        mw=pt.mw_c,
    )
