"""Bartz correlation: the equations of Huzel and Huang (1992), pp. 85-86."""

import unittest

from ceatofoam import heat


class TestBartz(unittest.TestCase):
    def test_prandtl_estimate(self):
        # Eq. (4-15) at gamma = 1.4
        self.assertAlmostEqual(heat.prandtl_estimate(1.4), 5.6/7.6)

    def test_recovery(self):
        self.assertAlmostEqual(heat.recovery_factor(0.81, turbulent=False),
                               0.9)
        self.assertAlmostEqual(heat.recovery_factor(0.7), 0.7**0.33)

    def test_adiabatic_wall_temperature(self):
        self.assertEqual(heat.adiabatic_wall_temperature(3000., 1.2, 0., .9),
                         3000.)
        # r = 1: the stagnation temperature at any Mach number
        self.assertAlmostEqual(
            heat.adiabatic_wall_temperature(3000., 1.2, 2.5, 1.0), 3000.)
        self.assertLess(
            heat.adiabatic_wall_temperature(3000., 1.2, 2.5, 0.9), 3000.)

    def test_sigma_is_one_at_the_limit(self):
        # T_w = T_c, M = 0: no property variation across the layer
        self.assertAlmostEqual(heat.sigma(3000., 3000., 1.2, 0.0), 1.0)
        # Fig. 4-28: a cold wall raises sigma above 1 at the throat
        self.assertGreater(heat.sigma(600., 3000., 1.2, 1.0), 1.0)

    def test_bartz_scales(self):
        g = heat.ChamberGas("test", 1e-4, 2000., 0.7, 1.2)
        h1 = heat.bartz(g, 3.5e6, 1750., 0.05, 0.0375, 1.0, 1.0)
        h2 = heat.bartz(g, 7.0e6, 1750., 0.05, 0.0375, 1.0, 1.0)
        h3 = heat.bartz(g, 3.5e6, 1750., 0.05, 0.0375, 4.0, 1.0)
        self.assertAlmostEqual(h2/h1, 2**0.8)
        self.assertAlmostEqual(h3/h1, 4**-0.9)
        expected = (0.026/0.05**0.2*(1e-4**0.2*2000./0.7**0.6)
                    *(3.5e6/1750.)**0.8*(0.05/0.0375)**0.1)
        self.assertAlmostEqual(h1, expected)


if __name__ == "__main__":
    unittest.main()
