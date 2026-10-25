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

"""Ambient conditions from the NASA Glenn Earth atmosphere model.

Source: NASA Glenn Research Center, "Earth Atmosphere Model, Metric Units",
https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html . The model has
three layers with boundaries at 11 000 m and 25 000 m (metres, not feet).
Temperatures are in deg C and pressures in kPa in the source; the source uses
273.1 in its own formulas, which is kept here. The returned temperature is
converted to kelvin with 273.15.
"""

import math


def ambient(altitude):
    """Return (temperature [K], pressure [Pa]) at ``altitude`` [m]."""
    h = float(altitude)
    if h < 11000.0:
        # Troposphere
        t = 15.04 - 0.00649*h
        p = 101.29*((t + 273.1)/288.08)**5.256
    elif h < 25000.0:
        # Lower stratosphere
        t = -56.46
        p = 22.65*math.exp(1.73 - 0.000157*h)
    else:
        # Upper stratosphere
        t = -131.21 + 0.00299*h
        p = 2.488*((t + 273.1)/216.6)**(-11.388)
    return t + 273.15, p*1.0e3
