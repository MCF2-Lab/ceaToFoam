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
