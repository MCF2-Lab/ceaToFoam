"""Conical thrust-chamber contour.

The wall is, from the injector face to the exit: a cylinder of radius rc, a
fillet arc of radius r_cf tangent to the cylinder and to the convergent cone
(half-angle theta_c), the convergent cone, an upstream throat arc of radius
r_tu, a downstream throat arc of radius r_td and the divergent cone
(half-angle theta_d) to the exit radius re. All arcs are tangent to their
neighbours. The cylinder length is chosen so that the volume from the
injector face to the throat is L* A_t (Huzel and Huang 1992, Sec. 4.3).

Angles are in degrees, lengths in metres.
"""

from dataclasses import dataclass, field
import math

import numpy as np


def contraction_ratio_correlation(throat_diameter):
    """Empirical chamber contraction ratio Ac/At = 8.0 Dt^-0.6 + 1.25.

    Dt in cm (Humble, Henry and Larson 1995, Space Propulsion Analysis and
    Design, Fig. 5.24, adapted from Huzel and Huang 1992).
    ``throat_diameter`` is given in metres.
    """
    return 8.0*(throat_diameter*100.0)**-0.6 + 1.25


@dataclass
class Segment:
    """One wall segment between two contour stations."""

    kind: str           # "line" or "arc"
    start: tuple
    end: tuple
    mid: tuple = None   # point on the arc halfway along it (arcs only)

    @property
    def length(self):
        if self.kind == "line":
            return math.dist(self.start, self.end)
        # Arc through start, mid, end: length = radius * angle.
        (x1, y1), (x2, y2), (x3, y3) = self.start, self.mid, self.end
        a, b, c = (
            math.dist(self.start, self.mid),
            math.dist(self.mid, self.end),
            math.dist(self.start, self.end),
        )
        area2 = abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
        radius = a*b*c/(2*area2)
        return 2*radius*math.asin(min(1.0, c/(2*radius)))


@dataclass
class Contour:
    """Thrust-chamber wall from the injector face (x = 0) to the exit."""

    rt: float
    rc: float
    re: float
    theta_c: float
    theta_d: float
    r_cf: float
    r_tu: float
    r_td: float
    cylinder_length: float
    segments: list = field(default_factory=list)

    @property
    def stations(self):
        """Wall points (x, r) at the segment ends, injector face first."""
        return [self.segments[0].start] + [s.end for s in self.segments]

    @property
    def x_throat(self):
        return self.segments[3].end[0]

    @property
    def x_exit(self):
        return self.segments[-1].end[0]

    @property
    def area_ratio(self):
        return (self.re/self.rt)**2

    @property
    def contraction_ratio(self):
        return (self.rc/self.rt)**2

    def radius(self, x):
        """Wall radius at axial position x (vectorized)."""
        x = np.asarray(x, dtype=float)
        r = np.full_like(x, np.nan)
        for seg in self.segments:
            (x0, y0), (x1, y1) = seg.start, seg.end
            if x1 <= x0:
                continue
            inside = (x >= x0) & (x <= x1)
            if seg.kind == "line":
                r[inside] = y0 + (y1 - y0)*(x[inside] - x0)/(x1 - x0)
            else:
                cx, cy, rad = _arc_centre(seg)
                sign = 1.0 if seg.mid[1] > cy else -1.0
                r[inside] = cy + sign*np.sqrt(
                    np.maximum(rad**2 - (x[inside] - cx)**2, 0.0)
                )
        return r

    def volume(self, x0, x1, n=4001):
        """Volume enclosed by the wall between x0 and x1."""
        x = np.linspace(x0, x1, n)
        return float(np.trapezoid(math.pi*self.radius(x)**2, x))


def _arc_centre(seg):
    """Centre and radius of the circle through an arc's three points."""
    (x1, y1), (x2, y2), (x3, y3) = seg.start, seg.mid, seg.end
    d = 2*(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))
    ux = (
        (x1**2 + y1**2)*(y2 - y3) + (x2**2 + y2**2)*(y3 - y1)
        + (x3**2 + y3**2)*(y1 - y2)
    )/d
    uy = (
        (x1**2 + y1**2)*(x3 - x2) + (x2**2 + y2**2)*(x1 - x3)
        + (x3**2 + y3**2)*(x2 - x1)
    )/d
    return ux, uy, math.dist((ux, uy), seg.start)


def _arc(centre, radius, a0, a1):
    """Wall arc on a circle; angles in radians measured from the centre."""
    cx, cy = centre
    pt = lambda a: (cx + radius*math.cos(a), cy + radius*math.sin(a))
    return Segment("arc", pt(a0), pt(a1), pt(0.5*(a0 + a1)))


def conical_contour(
    rt, rc, re, theta_c, theta_d, r_cf, r_tu, r_td, chamber_volume,
):
    """Build the wall and solve the cylinder length for ``chamber_volume``.

    ``chamber_volume`` is the volume from the injector face to the throat
    (L* A_t). Raises ValueError if the convergent section alone already
    exceeds it or if the radii are inconsistent.
    """
    if not rc > rt or not re > rt:
        raise ValueError("need rc > rt and re > rt")
    tc, td = math.radians(theta_c), math.radians(theta_d)

    # Throat at x = 0 for now; arcs below the centre (wall is convex
    # towards the axis at the throat, concave at the chamber fillet).
    # Upstream throat arc: centre (0, rt + r_tu), from the cone tangency
    # point at angle -(pi/2 + tc) to the throat at -pi/2.
    up = _arc((0.0, rt + r_tu), r_tu, -0.5*math.pi - tc, -0.5*math.pi)
    down = _arc((0.0, rt + r_td), r_td, -0.5*math.pi, -0.5*math.pi + td)
    xa, ya = up.start
    xb, yb = down.end
    if yb >= re:
        raise ValueError("downstream throat arc reaches the exit radius")

    # Chamber fillet: centre (xf, rc - r_cf); tangent to the cylinder at
    # (xf, rc) and to the cone at angle (pi/2 - tc) from the centre.
    y_cone = rc - r_cf*(1 - math.cos(tc))
    if y_cone <= ya:
        raise ValueError("chamber fillet too large for the convergent cone")
    x_cone = xa - (y_cone - ya)/math.tan(tc)
    xf = x_cone - r_cf*math.sin(tc)
    fillet = _arc((xf, rc - r_cf), r_cf, 0.5*math.pi, 0.5*math.pi - tc)
    cone_c = Segment("line", fillet.end, up.start)
    x_exit = xb + (re - yb)/math.tan(td)
    cone_d = Segment("line", down.end, (x_exit, re))

    trial = Contour(
        rt, rc, re, theta_c, theta_d, r_cf, r_tu, r_td, 0.0,
        [Segment("line", (xf, rc), (xf, rc)), fillet, cone_c, up, down,
         cone_d],
    )
    v_conv = trial.volume(xf, 0.0)
    v_cyl = chamber_volume - v_conv
    if v_cyl <= 0.0:
        raise ValueError(
            f"L* too small: the convergent section holds {v_conv:.4g} m3,"
            f" more than the chamber volume {chamber_volume:.4g} m3"
        )
    l_cyl = v_cyl/(math.pi*rc**2)

    shift = l_cyl - xf

    def move(seg):
        mv = lambda p: None if p is None else (p[0] + shift, p[1])
        return Segment(seg.kind, mv(seg.start), mv(seg.end), mv(seg.mid))

    segments = [Segment("line", (0.0, rc), (l_cyl, rc))] + [
        move(s) for s in (fillet, cone_c, up, down, cone_d)
    ]
    return Contour(
        rt, rc, re, theta_c, theta_d, r_cf, r_tu, r_td, l_cyl, segments,
    )
