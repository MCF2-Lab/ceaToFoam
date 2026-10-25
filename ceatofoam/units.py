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

"""Physical constants and parsing of dimensional inputs.

Inputs may be plain numbers (SI) or strings such as ``"500 psi"``,
``"2250 lbf"``, ``"40 in"`` or ``"15 deg"``. Everything inside the package is
SI, with angles in degrees.
"""

import re

#: Universal gas constant [J/(kmol K)], CODATA 2018 (exact).
RU = 8314.462618
#: Standard gravity [m/s^2] (exact, CGPM 1901).
G0 = 9.80665
#: Standard pressure for entropy [Pa] (NASA thermodynamic data, 1 bar).
P_STD = 1.0e5
#: Reference temperature for enthalpy [K].
T_REF = 298.15

# Conversion factors to SI (exact definitions).
_FACTORS = {
    # pressure -> Pa
    "pa": 1.0, "kpa": 1.0e3, "mpa": 1.0e6, "bar": 1.0e5,
    "atm": 101325.0, "psi": 6894.757293168361, "psia": 6894.757293168361,
    # force -> N
    "n": 1.0, "kn": 1.0e3, "mn": 1.0e6, "lbf": 4.4482216152605,
    # length -> m
    "m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3, "in": 0.0254, "ft": 0.3048,
    # temperature -> K
    "k": 1.0,
    # angle -> deg
    "deg": 1.0,
    # dimensionless
    "": 1.0,
}

_DIMENSION = {
    "pressure": {"pa", "kpa", "mpa", "bar", "atm", "psi", "psia"},
    "force": {"n", "kn", "mn", "lbf"},
    "length": {"m", "cm", "mm", "in", "ft"},
    "temperature": {"k"},
    "angle": {"deg"},
    "dimensionless": {""},
}

_PATTERN = re.compile(r"^\s*([-+0-9.eE]+)\s*([A-Za-z]*)\s*$")


def to_si(value, dimension):
    """Convert ``value`` (number or ``"<number> <unit>"``) to SI.

    ``dimension`` is one of pressure, force, length, temperature, angle or
    dimensionless; a number without a unit is SI, and a unit of another
    dimension raises ValueError.
    """
    if isinstance(value, (int, float)):
        return float(value)
    match = _PATTERN.match(str(value))
    if not match:
        raise ValueError(f"cannot parse quantity {value!r}")
    number, unit = float(match.group(1)), match.group(2).lower()
    if unit and unit not in _DIMENSION[dimension]:
        raise ValueError(f"unit {unit!r} in {value!r} is not a {dimension}")
    return number*_FACTORS[unit]
