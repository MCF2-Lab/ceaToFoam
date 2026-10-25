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

"""Units, atmosphere and constant-gamma relations (no external programs)."""

import math
import unittest

from ceatofoam import atmosphere, sizing
from ceatofoam.units import to_si


class TestUnits(unittest.TestCase):
    def test_conversions(self):
        self.assertAlmostEqual(to_si("500 psi", "pressure"), 3447378.6466,
                               places=3)
        self.assertAlmostEqual(to_si("1 lbf", "force"), 4.4482216152605)
        self.assertAlmostEqual(to_si("40 in", "length"), 1.016)
        self.assertEqual(to_si(15, "angle"), 15.0)

    def test_number_string_is_si(self):
        # the CLI passes strings, e.g. the default --altitude "0"
        self.assertEqual(to_si("0", "length"), 0.0)
        self.assertEqual(to_si("3.5e6", "pressure"), 3.5e6)

    def test_wrong_dimension(self):
        with self.assertRaises(ValueError):
            to_si("500 psi", "length")


class TestAtmosphere(unittest.TestCase):
    def test_sea_level(self):
        T, p = atmosphere.ambient(0.0)
        self.assertAlmostEqual(T, 288.19, places=2)
        self.assertAlmostEqual(p, 101.29e3*(288.14/288.08)**5.256, places=3)

    def test_layer_boundaries_in_metres(self):
        # ceaToFoam 1.x placed the boundaries at 11 000 ft and 25 000 ft.
        # At 5000 m the troposphere formula applies.
        T, p = atmosphere.ambient(5000.0)
        self.assertAlmostEqual(T, 15.04 - 0.00649*5000 + 273.15, places=9)
        self.assertAlmostEqual(
            p, 101.29e3*((15.04 - 0.00649*5000 + 273.1)/288.08)**5.256
        )

    def test_continuity(self):
        for h in (11000.0, 25000.0):
            _, below = atmosphere.ambient(h - 1e-6)
            _, above = atmosphere.ambient(h)
            self.assertLess(abs(above/below - 1), 0.02)


class TestConstantGamma(unittest.TestCase):
    def test_area_ratio(self):
        # A/A* = 1.6875 at M = 2, gamma = 1.4 (NACA 1135).
        self.assertAlmostEqual(
            sizing.area_ratio_from_mach(2.0, 1.4), 1.6875, places=4
        )

    def test_mass_flux_function(self):
        self.assertAlmostEqual(sizing.mass_flux_function(1.4), 0.684731,
                               places=6)

    def test_reproduces_v1_appendix(self):
        # ceaToFoam 1.x sizing (paper Appendix A): gamma_s = 1.1389,
        # averaged M = 23.12, Tc = 3526.5 K, 500 psi, 500 lbf, sea level.
        d = sizing.size_constant_gamma(
            thrust=500*4.4482216152605, pc=500*6894.757293168361,
            tc=3526.5, p_amb=101400.93090454885, gamma=1.1389, mw=23.12,
        )
        self.assertAlmostEqual(d.mdot/0.8249168314057205, 1, places=5)
        self.assertAlmostEqual(d.area_ratio/5.878210176427288, 1, places=5)
        self.assertAlmostEqual(d.at/0.0004234403665069037, 1, places=5)
        self.assertAlmostEqual(d.mach_e/2.781606421750013, 1, places=5)


if __name__ == "__main__":
    unittest.main()
