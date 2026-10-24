"""Frozen mixture, janaf and Sutherland fits, frozen sizing against CEA."""

import os
import shutil
import unittest

import numpy as np

try:
    import cantera  # noqa: F401
    HAVE_CANTERA = True
except ImportError:
    HAVE_CANTERA = False


def have_cea():
    exe = os.environ.get("CEA_EXE_LOCATION")
    if exe:
        return os.access(exe, os.X_OK)
    try:
        import CEA_Wrap  # noqa: F401
    except ImportError:
        return False
    return os.name == "nt" or shutil.which("Dan_CEA2") is not None


AIR = {"N2": 0.78084, "O2": 0.20946, "AR": 0.00934}


@unittest.skipUnless(HAVE_CANTERA, "Cantera not installed")
class TestMixture(unittest.TestCase):
    def setUp(self):
        from ceatofoam.thermo import FrozenMixture

        self.air = FrozenMixture(AIR)

    def test_air(self):
        self.assertAlmostEqual(self.air.mw, 28.96, places=2)
        self.assertAlmostEqual(self.air.gamma(300.0), 1.400, delta=1e-3)
        self.assertEqual(self.air.dropped_fraction, 0.0)

    def test_dropped_species(self):
        from ceatofoam.thermo import FrozenMixture

        m = FrozenMixture({"N2": 0.99, "XYZ": 0.01})
        self.assertAlmostEqual(m.dropped_fraction, 0.01)
        self.assertAlmostEqual(m.composition["N2"], 1.0)

    def test_janaf_continuity_and_error(self):
        from ceatofoam.thermo import JanafFit, fit_janaf

        fit = fit_janaf(self.air)
        tc = fit.t_common
        for f in (JanafFit.cp_over_r, JanafFit.h_over_r, JanafFit.s_over_r):
            lo, hi = f(fit.low, tc), f(fit.high, tc)
            self.assertLess(abs(hi - lo), 1e-9*abs(lo))
        self.assertLess(fit.cp_max_rel_error, 2e-3)
        # Enthalpy at T_ref reproduces the mixture value.
        self.assertAlmostEqual(
            JanafFit.h_over_r(fit.low, 298.15)*self.air.R,
            self.air.h(298.15), delta=1.0,
        )

    def test_sutherland(self):
        from ceatofoam.thermo import fit_sutherland

        s = fit_sutherland(self.air, 250.0, 1500.0)
        self.assertLess(s.max_rel_error, 0.03)
        # Air: mu(300 K) = 1.85e-5 Pa s within the fit error.
        self.assertAlmostEqual(s.mu(300.0)/1.85e-5, 1, delta=0.03)

    def test_frozen_expansion_air(self):
        """Frozen expansion of nearly calorically perfect cold air agrees
        with the constant-gamma relations."""
        from ceatofoam import sizing

        tc, pc = 300.0, 5.0e5
        st = sizing.frozen_throat(self.air, tc, pc)
        g = self.air.gamma(0.5*(tc + st.T))
        self.assertAlmostEqual(st.p/pc, (2/(g + 1))**(g/(g - 1)), places=3)


@unittest.skipUnless(HAVE_CANTERA and have_cea(), "CEA executable missing")
class TestAgainstCEA(unittest.TestCase):
    """The frozen sizing reproduces CEA's own frozen (nfz = 1) solution."""

    @classmethod
    def setUpClass(cls):
        from ceatofoam import cea, thermo

        cls.pc = 500*6894.757293168361
        cls.p_amb = 101400.93090454885
        cls.frozen = cea.run_rocket(
            "RP-1", "O2(L)", cls.pc, 2.4, cls.pc/cls.p_amb, "frozen"
        )
        cls.mix = thermo.FrozenMixture(cls.frozen.composition)

    def test_chamber_properties(self):
        f = self.frozen
        self.assertAlmostEqual(self.mix.gamma(f.tc)/f.gamma_c, 1, places=3)
        self.assertAlmostEqual(self.mix.cp(f.tc)/f.cp_c, 1, places=3)
        self.assertAlmostEqual(self.mix.mw/f.mw_c, 1, places=4)

    def test_v1_station_average_molecular_weight(self):
        """ceaToFoam 1.x averaged M,(1/n) over chamber, throat and exit:
        23.12 kg/kmol for this case (paper Appendix A)."""
        from ceatofoam import cea

        eq = cea.run_rocket(
            "RP-1", "O2(L)", self.pc, 2.4, self.pc/self.p_amb, "equilibrium"
        )
        self.assertAlmostEqual(sum(eq.m_stations)/3, 23.12, places=2)

    def test_exit_state(self):
        from ceatofoam import sizing

        f = self.frozen
        d = sizing.size_frozen(2224.11, self.pc, f.tc, self.p_amb, self.mix)
        self.assertAlmostEqual(d.te/f.t_e, 1, places=3)
        self.assertAlmostEqual(d.mach_e/f.mach_e, 1, places=3)
        self.assertAlmostEqual(d.isp/f.isp, 1, places=3)
        self.assertAlmostEqual(d.cstar/f.cstar, 1, places=3)


if __name__ == "__main__":
    unittest.main()
