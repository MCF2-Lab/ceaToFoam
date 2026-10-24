"""ceaToFoam: NASA CEA chamber thermochemistry to OpenFOAM rocket-nozzle cases.

The package sizes a conical thrust chamber from a CEA chamber solution and
writes an axisymmetric rhoCentralFoam case whose thermophysical model is the
frozen chamber mixture. Sizing and CFD use the same frozen thermodynamics
unless another sizing model is requested explicitly.
"""

__version__ = "2.0.0.dev0"
