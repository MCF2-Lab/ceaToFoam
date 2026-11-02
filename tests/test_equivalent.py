"""Equivalent gas: the first law on the path, and the expansion against CEA."""

import math
import os
import shutil
import unittest

import numpy as np

from ceatofoam import equivalent, sizing, thermo
from ceatofoam.geometry import conical_contour

try:
    import cantera  # noqa: F401
    HAVE_CANTERA = True
except ImportError:
    HAVE_CANTERA = False


def have_cea():
    exe = os.environ.get("CEA_EXE_LOCATION")
    if exe:
        return os.access(exe, os.X_OK)
    return shutil.which("Dan_CEA2") is not None


# Products of RP-1/LOX at 500 psi, O/F 2.4 (rounded): any close state serves,
# the tests equilibrate it at the chamber temperature.
COMPOSITION = {"CO": 0.33, "CO2": 0.11, "H2": 0.12, "H2O": 0.32, "H": 0.02,
               "O": 0.01, "OH": 0.07, "O2": 0.02}
TC, PC, PA = 3526.5, 3.447e6, 1.014e5


@unittest.skipUnless(HAVE_CANTERA, "needs Cantera")
class TestPath(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mix = thermo.FrozenMixture(COMPOSITION)
        cls.path = equivalent.equilibrium_path(cls.mix.composition, TC, PC,
                                               0.15*PA)
        cls.gas = equivalent.EquivalentGas(cls.path, cls.mix, TC)

    def test_chamber_state(self):
        # tau equals T at the chamber, by the choice of the gas constant
        self.assertAlmostEqual(self.gas.tau[-1]/TC, 1.0, places=9)

    def test_first_law_on_the_path(self):
        p, rho, e = self.path.p, self.path.rho, self.path.e
        lhs = np.gradient(e, rho)[3:-3]
        rhs = (p/rho**2)[3:-3]
        self.assertLess(np.max(np.abs(lhs/rhs - 1)), 2e-4)

    def test_sound_speed_of_the_gas_is_that_of_the_path(self):
        p, rho = self.path.p, self.path.rho
        a2 = np.gradient(p, rho)[5:-5]
        tau = (p/(rho*self.gas.R))[5:-5]
        cv = np.interp(tau, self.gas.tau, self.gas.cv_path)
        a2_gas = (cv + self.gas.R)/cv*self.gas.R*tau
        self.assertLess(np.max(np.abs(a2_gas/a2 - 1)), 1e-3)

    def test_fit(self):
        self.assertLess(self.gas.janaf.cp_max_rel_error, 5e-3)
        t_lo = self.gas.janaf.t_low
        t_hi = self.gas.janaf.t_high
        for t in np.linspace(t_lo, t_hi, 200):
            self.assertGreater(self.gas.cp(t)/self.gas.R, 1.0)

    def test_expansion_of_the_fitted_gas_is_the_path(self):
        # The frozen sizing functions, run on the fitted gas, against the
        # path itself (no CEA here).
        one_d = sizing.size_frozen(1.0, PC, TC, PA, self.gas)
        u = self.path.velocity
        ie = int(np.argmin(np.abs(self.path.p - PA)))
        g = self.path.mass_flux
        self.assertAlmostEqual(one_d.ve/u[ie], 1.0, delta=1e-3)
        self.assertAlmostEqual(one_d.area_ratio/(g.max()/g[ie]), 1.0,
                               delta=2e-3)
        self.assertAlmostEqual(one_d.cstar/(PC/g.max()), 1.0, delta=1e-3)

    def test_true_temperature(self):
        self.assertAlmostEqual(float(self.gas.true_temperature(TC)), TC,
                               delta=0.5)
        ie = int(np.argmin(np.abs(self.path.p - PA)))
        tau_e = self.path.p[ie]/(self.path.rho[ie]*self.gas.R)
        self.assertAlmostEqual(float(self.gas.true_temperature(tau_e)),
                               self.path.T[ie], delta=1.0)


@unittest.skipUnless(HAVE_CANTERA, "needs Cantera")
class TestKinetic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mix = thermo.FrozenMixture(COMPOSITION)
        cls.path_eq = equivalent.equilibrium_path(
            cls.mix.composition, TC, PC, 0.15*PA)
        cls.uf = cls._frozen_exit_velocity()

    @classmethod
    def _frozen_exit_velocity(cls):
        return sizing.frozen_station(cls._frozen(), TC, PC, PA).u

    @classmethod
    def _frozen(cls):
        import cantera as ct
        gas = ct.Solution("gri30_highT.yaml")
        gas.TPX = TC, PC, cls.mix.composition
        gas.equilibrate("TP")
        comp = {s: x for s, x in zip(gas.species_names, gas.X) if x > 1e-9}
        return thermo.FrozenMixture(comp)

    def _exit_velocity(self, rt):
        rc, re = rt*2.5, rt*math.sqrt(5.0)
        contour = conical_contour(rt, rc, re, 30.0, 15.0, rt, 1.5*rt, 1.5*rt,
                                  40*math.pi*rt**3)
        path = equivalent.kinetic_path(self.mix.composition, TC, PC, 0.15*PA,
                                       contour, n=200)
        ie = int(np.argmin(np.abs(path.p - PA)))
        return float(path.velocity[ie])

    def test_between_the_bounds_and_rising_with_scale(self):
        ie = int(np.argmin(np.abs(self.path_eq.p - PA)))
        u_eq = float(self.path_eq.velocity[ie])
        small, large = self._exit_velocity(0.0025), self._exit_velocity(0.25)
        self.assertGreater(small, self.uf)
        self.assertLess(small, large)
        self.assertLess(large, u_eq*(1 + 1e-4))


if __name__ == "__main__":
    unittest.main()
