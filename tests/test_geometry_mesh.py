"""Contour construction and the wedge blockMeshDict."""

import math
import os
import re
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

from ceatofoam import geometry, mesh


def make_contour(theta_c=30.0, theta_d=15.0):
    rt = 0.0116
    at = math.pi*rt**2
    return geometry.conical_contour(
        rt, rt*math.sqrt(6.0), rt*math.sqrt(5.0), theta_c, theta_d,
        rt, 1.5*rt, 1.5*rt, 1.016*at,
    )


def slope(seg, at_start):
    """dr/dx of a segment at its start or end."""
    if seg.kind == "line":
        (x0, y0), (x1, y1) = seg.start, seg.end
        return (y1 - y0)/(x1 - x0)
    cx, cy, _ = geometry._arc_centre(seg)
    x, y = seg.start if at_start else seg.end
    return -(x - cx)/(y - cy)


class TestContour(unittest.TestCase):
    def setUp(self):
        self.c = make_contour()

    def test_continuity_and_tangency(self):
        segs = self.c.segments
        for a, b in zip(segs[:-1], segs[1:]):
            self.assertLess(math.dist(a.end, b.start), 1e-12)
            self.assertAlmostEqual(slope(a, False), slope(b, True), places=9)

    def test_angles_and_radii(self):
        s = self.c.segments
        self.assertAlmostEqual(slope(s[2], True),
                               -math.tan(math.radians(30.0)), places=12)
        self.assertAlmostEqual(slope(s[5], True),
                               math.tan(math.radians(15.0)), places=12)
        self.assertAlmostEqual(s[3].end[1], self.c.rt, places=15)
        self.assertAlmostEqual(s[5].end[1], self.c.re, places=15)
        x = np.linspace(0, self.c.x_exit, 20001)
        self.assertAlmostEqual(np.nanmin(self.c.radius(x)), self.c.rt,
                               places=8)

    def test_chamber_volume(self):
        at = math.pi*self.c.rt**2
        v = self.c.volume(0.0, self.c.x_throat, n=200001)
        self.assertAlmostEqual(v/(1.016*at), 1, places=6)

    def test_lstar_too_small(self):
        rt = 0.0116
        with self.assertRaises(ValueError):
            geometry.conical_contour(
                rt, 3*rt, 2*rt, 30, 15, rt, 1.5*rt, 1.5*rt, 1e-9
            )


class TestWedge(unittest.TestCase):
    def setUp(self):
        self.c = make_contour()
        self.text, self.summary = mesh.block_mesh_dict(
            self.c, mesh.MeshSettings()
        )
        body = self.text.split("vertices")[1].split(");")[0]
        self.verts = np.array([
            [float(v) for v in m.groups()]
            for m in re.finditer(r"\(([-0-9.e+]+) ([-0-9.e+]+) ([-0-9.e+]+)\)",
                                 body)
        ])

    def test_symmetric_and_radius_preserved(self):
        v = self.verts
        off = v[:, 1] > 0
        r = np.hypot(v[off, 1], v[off, 2])
        angle = np.degrees(np.arctan2(np.abs(v[off, 2]), v[off, 1]))
        self.assertTrue(np.allclose(angle, 2.5, atol=1e-12))
        # Every off-axis vertex has a mirror image.
        plus = {(round(x, 12), round(y, 12), round(z, 12)) for x, y, z in v}
        for x, y, z in v[off]:
            self.assertIn((round(x, 12), round(y, 12), round(-z, 12)), plus)
        # Wall vertices keep the contour radius.
        for x, rr in self.c.stations:
            hit = np.isclose(v[:, 0], x) & np.isclose(r_full(v), rr)
            self.assertTrue(hit.any())
        self.assertTrue(np.all(r > 0))

    def test_resolution_scales(self):
        _, s2 = mesh.block_mesh_dict(self.c, mesh.MeshSettings(scale=2.0))
        self.assertEqual(s2["radial_cells"], 2*self.summary["radial_cells"])
        # Chamber cells scale with scale^2 (up to rounding of the axial
        # counts); the geometric plume cells grow only logarithmically.
        nozzle = lambda s: s["radial_cells"]*sum(s["axial_cells"])
        ratio = nozzle(s2)/nozzle(self.summary)
        self.assertGreater(ratio, 3.9)
        self.assertLess(ratio, 4.1)

    @unittest.skipUnless(shutil.which("blockMesh") and shutil.which("checkMesh"),
                         "OpenFOAM not sourced")
    def test_check_mesh(self):
        tut = os.path.join(os.environ["FOAM_TUTORIALS"], "compressible",
                           "rhoCentralFoam", "LadenburgJet60psi", "system")
        with tempfile.TemporaryDirectory() as case:
            os.makedirs(os.path.join(case, "system"))
            for f in ("controlDict", "fvSchemes", "fvSolution"):
                shutil.copy(os.path.join(tut, f),
                            os.path.join(case, "system", f))
            mesh.write_block_mesh_dict(case, self.c, mesh.MeshSettings())
            for app in ("blockMesh", "checkMesh"):
                out = subprocess.run([app, "-case", case],
                                     capture_output=True, text=True)
                self.assertEqual(out.returncode, 0, out.stdout[-2000:])
            self.assertIn("Mesh OK.", out.stdout)


def r_full(v):
    return np.hypot(v[:, 1], v[:, 2])


if __name__ == "__main__":
    unittest.main()
