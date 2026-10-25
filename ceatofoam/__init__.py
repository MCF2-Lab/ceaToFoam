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

"""ceaToFoam: NASA CEA chamber thermochemistry to OpenFOAM rocket-nozzle cases.

The package sizes a conical thrust chamber from a CEA chamber solution and
writes an axisymmetric rhoCentralFoam case whose thermophysical model is the
frozen chamber mixture. Sizing and CFD use the same frozen thermodynamics
unless another sizing model is requested explicitly.
"""

__version__ = "2.0.0"
