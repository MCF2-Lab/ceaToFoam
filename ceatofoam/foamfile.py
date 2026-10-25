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

"""Minimal writer for OpenFOAM dictionary files."""

import os

_BANNER = r"""/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Written by ceaToFoam
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
"""

_SEPARATOR = (
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
    "* * * //\n"
)
_END = (
    "\n// *************************************************************"
    "************ //\n"
)


def header(cls, obj, location=None):
    lines = [
        "FoamFile\n{\n",
        "    version     2.0;\n",
        "    format      ascii;\n",
        f"    class       {cls};\n",
    ]
    if location:
        lines.append(f'    location    "{location}";\n')
    lines.append(f"    object      {obj};\n}}\n")
    return _BANNER + "".join(lines) + _SEPARATOR + "\n"


def write(case, relpath, cls, body):
    """Write ``body`` with a FoamFile header to ``case/relpath``."""
    path = os.path.join(case, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    location, obj = os.path.split(relpath)
    with open(path, "w") as f:
        f.write(header(cls, obj, location or None))
        f.write(body.rstrip() + "\n")
        f.write(_END)
    return path


def fmt(x):
    """Shortest representation that round-trips the double exactly."""
    return repr(float(x))
