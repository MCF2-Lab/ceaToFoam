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

"""Axisymmetric wedge blockMeshDict for the thrust chamber and plume.

The wedge has half-angle alpha on either side of the x-y plane: a point at
axial position x and radius r gives the vertices (x, r cos(alpha),
+-r sin(alpha)), so both wedge patches are mirror images and every vertex
keeps its radius. Each wall segment of the contour is one block between the
axis and the wall; two further blocks hold the plume behind the exit plane.

Resolution: ``radial_cells`` cells between the axis and the wall at every
axial station, and axial cells of about ``cell_aspect`` times the local
radial cell size, so the cell count follows the local radius and the throat
is always resolved. ``scale`` multiplies every cell count (mesh studies). The
cylindrical chamber is meshed only with ``include_cylinder``: in a steady
frozen flow it is a duct at constant total conditions, and its low Mach
number made the density-based solver converge slowly with spurious
oscillations near the axis; by default the inlet is the end of the
cylinder.
By default the domain ends at the exit plane (supersonic outflow), which is
all that thrust and mass flow need. With ``plume_length`` > 0 a plume region
extends that many exit diameters downstream and ``plume_radius`` exit radii
radially; its cells grow geometrically from the exit-plane size with ratio
``growth``.

Patches: inlet (injector face), nozzle (wall), outlet (the exit plane, or
the downstream plane of the plume region), front and back (wedge). With a
plume region also farfield (outer radius) and base (wall: the exit-plane
annulus outside the lip, the engine base; an open boundary there touches
the lip corner, where the expanding jet drove the temperature out of range
in tests).
"""

from dataclasses import dataclass
import math

from .foamfile import fmt, write


@dataclass
class MeshSettings:
    radial_cells: int = 40
    cell_aspect: float = 1.0
    scale: float = 1.0
    wedge_angle: float = 5.0      # total angle [deg]
    plume_length: float = 0.0     # exit diameters; 0: no plume region
    plume_radius: float = 3.0     # exit radii
    growth: float = 1.05
    include_cylinder: bool = False   # mesh the cylindrical chamber too


def _geometric_cells(first, length, growth):
    """Cell count and blockMesh expansion ratio (last/first) for cells that
    start at size ``first`` and grow by ``growth`` over ``length``."""
    if growth <= 1.0 + 1e-12:
        n = max(1, round(length/first))
        return n, 1.0
    n = math.log(1 + length*(growth - 1)/first)/math.log(growth)
    n = max(1, round(n))
    return n, growth**(n - 1)


def block_mesh_dict(contour, settings):
    """Return (text of system/blockMeshDict, summary dict)."""
    s = settings
    alpha = math.radians(0.5*s.wedge_angle)
    ca, sa = math.cos(alpha), math.sin(alpha)
    nr = max(4, round(s.radial_cells*s.scale))

    stations = contour.stations           # 7 wall points (x, r)
    x_e, r_e = stations[-1]
    plume = s.plume_length > 0.0
    x_end = x_e + s.plume_length*2*r_e if plume else x_e
    r_far = s.plume_radius*r_e if plume else r_e

    verts = []

    def add(x, r, side):
        verts.append((x, r*ca, side*r*sa))
        return len(verts) - 1

    axis = [add(x, 0.0, 0) for x, _ in stations]
    wm = [add(x, r, -1) for x, r in stations]
    wp = [add(x, r, +1) for x, r in stations]

    blocks, edges = [], []
    n_axial = []
    first = 0 if s.include_cylinder else 1
    for i, seg in enumerate(contour.segments):
        if i < first:
            n_axial.append(0)
            continue
        r_mean = 0.5*(seg.start[1] + seg.end[1])
        dx = s.cell_aspect*r_mean/nr
        n = max(2, math.ceil(seg.length/dx))
        n_axial.append(n)
        blocks.append(
            f"    hex ({axis[i]} {wm[i]} {wp[i]} {axis[i]} "
            f"{axis[i+1]} {wm[i+1]} {wp[i+1]} {axis[i+1]}) "
            f"({nr} 1 {n}) simpleGrading (1 1 1)"
        )
        if seg.kind == "arc":
            xm, rm = seg.mid
            for side, ids in ((-1, wm), (+1, wp)):
                edges.append(
                    f"    arc {ids[i]} {ids[i+1]} "
                    f"({fmt(xm)} {fmt(rm*ca)} {fmt(side*rm*sa)})"
                )

    used = range(first, 6)
    front = [f"({axis[i]} {wp[i]} {wp[i+1]} {axis[i+1]})" for i in used]
    back = [f"({axis[i]} {axis[i+1]} {wm[i+1]} {wm[i]})" for i in used]
    nozzle = [f"({wm[i]} {wp[i]} {wp[i+1]} {wm[i+1]})" for i in used]
    patches = [
        ("inlet", "patch",
         [f"({axis[first]} {wp[first]} {wm[first]} {axis[first]})"]),
        ("nozzle", "wall", nozzle),
    ]
    n_p = n_o = 0
    if plume:
        a7 = add(x_end, 0.0, 0)
        em, ep = add(x_end, r_e, -1), add(x_end, r_e, +1)
        gm, gp = add(x_e, r_far, -1), add(x_e, r_far, +1)
        fm, fp = add(x_end, r_far, -1), add(x_end, r_far, +1)
        # First axial cell equal to the last nozzle cell; first radial cell
        # outside the lip equal to the radial cell at the exit.
        dx_e = contour.segments[-1].length/n_axial[-1]
        n_p, g_p = _geometric_cells(dx_e, x_end - x_e, s.growth)
        n_o, g_o = _geometric_cells(r_e/nr, r_far - r_e, s.growth)
        blocks.append(
            f"    hex ({axis[6]} {wm[6]} {wp[6]} {axis[6]} "
            f"{a7} {em} {ep} {a7}) "
            f"({nr} 1 {n_p}) simpleGrading (1 1 {fmt(g_p)})"
        )
        blocks.append(
            f"    hex ({wm[6]} {gm} {gp} {wp[6]} {em} {fm} {fp} {ep}) "
            f"({n_o} 1 {n_p}) simpleGrading ({fmt(g_o)} 1 {fmt(g_p)})"
        )
        front += [f"({axis[6]} {wp[6]} {ep} {a7})", f"({wp[6]} {gp} {fp} {ep})"]
        back += [f"({axis[6]} {a7} {em} {wm[6]})", f"({wm[6]} {em} {fm} {gm})"]
        patches += [
            ("outlet", "patch",
             [f"({a7} {em} {ep} {a7})", f"({em} {fm} {fp} {ep})"]),
            ("farfield", "patch", [f"({gm} {gp} {fp} {fm})"]),
            ("base", "wall", [f"({wm[6]} {wp[6]} {gp} {gm})"]),
        ]
    else:
        patches.append(
            ("outlet", "patch", [f"({axis[6]} {wm[6]} {wp[6]} {axis[6]})"])
        )
    patches += [("front", "wedge", front), ("back", "wedge", back)]

    vtext = "\n".join(
        f"    ({fmt(x)} {fmt(y)} {fmt(z)}) // {k}"
        for k, (x, y, z) in enumerate(verts)
    )
    btext = "\n".join(
        f"    {name}\n    {{\n        type {kind};\n        faces\n        (\n"
        + "\n".join(f"            {f}" for f in faces)
        + "\n        );\n    }"
        for name, kind, faces in patches
    )
    text = f"""scale   1;

vertices
(
{vtext}
);

blocks
(
{chr(10).join(blocks)}
);

edges
(
{chr(10).join(edges)}
);

boundary
(
{btext}
);
"""
    n_cells = nr*sum(n_axial) + (nr + n_o)*n_p
    summary = {
        "radial_cells": nr,
        "axial_cells": n_axial,
        "plume": plume,
        "plume_axial_cells": n_p,
        "plume_radial_cells": n_o,
        "cells": n_cells,
        "x_exit": x_e,
        "x_inlet": stations[first][0],
        "x_end": x_end,
        "r_far": r_far,
        "wedge_angle": s.wedge_angle,
    }
    return text, summary


def write_block_mesh_dict(case, contour, settings):
    text, summary = block_mesh_dict(contour, settings)
    write(case, "system/blockMeshDict", "dictionary", text)
    return summary
