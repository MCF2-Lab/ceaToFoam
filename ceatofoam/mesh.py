"""Axisymmetric wedge blockMeshDict for the thrust chamber and plume.

The wedge has half-angle alpha on either side of the x-y plane: a point at
axial position x and radius r gives the vertices (x, r cos(alpha),
+-r sin(alpha)), so both wedge patches are mirror images and every vertex
keeps its radius. Each wall segment of the contour is one block between the
axis and the wall; two further blocks hold the plume behind the exit plane.

Resolution: ``radial_cells`` cells between the axis and the wall at every
axial station, and axial cells of about ``cell_aspect`` times the local
radial cell size, so the cell count follows the local radius and the throat
is always resolved. ``scale`` multiplies every cell count (mesh studies).
The plume extends ``plume_length`` exit diameters downstream and
``plume_radius`` exit radii radially; its cells grow geometrically from the
exit-plane size with ratio ``growth``.

Patches: inlet (injector face), nozzle (wall), base (wall: the exit-plane
annulus outside the lip, i.e. the engine base), outlet (downstream plane),
farfield (outer radius), front and back (wedge). An open boundary at the
base touches the lip corner, where the expanding jet drove the temperature
out of range in tests; the base is therefore a wall.
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
    plume_length: float = 4.0     # exit diameters
    plume_radius: float = 3.0     # exit radii
    growth: float = 1.05


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
    x_end = x_e + s.plume_length*2*r_e
    r_far = s.plume_radius*r_e

    verts = []

    def add(x, r, side):
        verts.append((x, r*ca, side*r*sa))
        return len(verts) - 1

    axis = [add(x, 0.0, 0) for x, _ in stations] + [add(x_end, 0.0, 0)]
    wm = [add(x, r, -1) for x, r in stations]
    wp = [add(x, r, +1) for x, r in stations]
    em, ep = add(x_end, r_e, -1), add(x_end, r_e, +1)
    gm, gp = add(x_e, r_far, -1), add(x_e, r_far, +1)
    fm, fp = add(x_end, r_far, -1), add(x_end, r_far, +1)

    blocks, edges = [], []
    n_axial = []
    for i, seg in enumerate(contour.segments):
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

    # Plume: first axial cell equal to the last nozzle cell; first radial
    # cell outside the lip equal to the radial cell at the exit.
    seg_e = contour.segments[-1]
    dx_e = seg_e.length/n_axial[-1]
    n_p, g_p = _geometric_cells(dx_e, x_end - x_e, s.growth)
    n_o, g_o = _geometric_cells(r_e/nr, r_far - r_e, s.growth)
    blocks.append(
        f"    hex ({axis[6]} {wm[6]} {wp[6]} {axis[6]} "
        f"{axis[7]} {em} {ep} {axis[7]}) "
        f"({nr} 1 {n_p}) simpleGrading (1 1 {fmt(g_p)})"
    )
    blocks.append(
        f"    hex ({wm[6]} {gm} {gp} {wp[6]} {em} {fm} {fp} {ep}) "
        f"({n_o} 1 {n_p}) simpleGrading ({fmt(g_o)} 1 {fmt(g_p)})"
    )

    front = [f"({axis[i]} {wp[i]} {wp[i+1]} {axis[i+1]})" for i in range(6)]
    back = [f"({axis[i]} {axis[i+1]} {wm[i+1]} {wm[i]})" for i in range(6)]
    front += [f"({axis[6]} {wp[6]} {ep} {axis[7]})", f"({wp[6]} {gp} {fp} {ep})"]
    back += [f"({axis[6]} {axis[7]} {em} {wm[6]})", f"({wm[6]} {em} {fm} {gm})"]
    nozzle = [f"({wm[i]} {wp[i]} {wp[i+1]} {wm[i+1]})" for i in range(6)]

    def faces(lst):
        return "\n".join(f"            {f}" for f in lst)

    vtext = "\n".join(
        f"    ({fmt(x)} {fmt(y)} {fmt(z)}) // {k}"
        for k, (x, y, z) in enumerate(verts)
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
    inlet
    {{
        type patch;
        faces
        (
            ({axis[0]} {wp[0]} {wm[0]} {axis[0]})
        );
    }}
    nozzle
    {{
        type wall;
        faces
        (
{faces(nozzle)}
        );
    }}
    outlet
    {{
        type patch;
        faces
        (
            ({axis[7]} {em} {ep} {axis[7]})
            ({em} {fm} {fp} {ep})
        );
    }}
    farfield
    {{
        type patch;
        faces
        (
            ({gm} {gp} {fp} {fm})
        );
    }}
    base
    {{
        type wall;
        faces
        (
            ({wm[6]} {wp[6]} {gp} {gm})
        );
    }}
    front
    {{
        type wedge;
        faces
        (
{faces(front)}
        );
    }}
    back
    {{
        type wedge;
        faces
        (
{faces(back)}
        );
    }}
);
"""
    n_cells = nr*sum(n_axial) + (nr + n_o)*n_p
    summary = {
        "radial_cells": nr,
        "axial_cells": n_axial,
        "plume_axial_cells": n_p,
        "plume_radial_cells": n_o,
        "cells": n_cells,
        "x_exit": x_e,
        "x_end": x_end,
        "r_far": r_far,
        "wedge_angle": s.wedge_angle,
    }
    return text, summary


def write_block_mesh_dict(case, contour, settings):
    text, summary = block_mesh_dict(contour, settings)
    write(case, "system/blockMeshDict", "dictionary", text)
    return summary
