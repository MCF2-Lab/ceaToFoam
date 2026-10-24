# ceaToFoam

ceaToFoam turns a NASA CEA chamber solution into a ready-to-run OpenFOAM
case for the compressible flow through a conical liquid-rocket thrust
chamber and its near plume. It

1. runs CEA (through [CEA_Wrap](https://github.com/civilwargeeky/CEA_Wrap))
   for the propellants, chamber pressure and mixture ratio, with shifting
   equilibrium and frozen at the chamber;
2. freezes the chamber composition and evaluates it with
   [Cantera](https://cantera.org) (molecular weight, cp(T), viscosity);
3. sizes the throat and exit for a design thrust and ambient pressure;
4. builds a conical chamber contour (cylinder, fillet, convergent cone,
   throat arcs, divergent cone) with the chamber volume given by L*;
5. writes an axisymmetric wedge `blockMeshDict` and a `rhoCentralFoam`
   case whose `janaf` and `sutherland` coefficients are fitted to the same
   frozen mixture;
6. after the run, evaluates thrust, mass flow and the exit state.

The sizing and the CFD use the same thermodynamics: by default the nozzle
is sized by the isentropic expansion of the frozen mixture with
temperature-dependent cp, which reproduces CEA's own frozen solution
(exit temperature, Mach number and Isp to within 0.1 %). Sizing with CEA's
shifting-equilibrium performance or with the constant-gamma relations is
available for comparison.

## Installation

Python 3.9 or newer and an OpenFOAM.com release (tested with v2606):

```sh
git clone https://github.com/MCF2-Lab/ceaToFoam.git
cd ceaToFoam
pip install .            # or: pip install .[plot] for the O/F sweep figure
```

CEA_Wrap ships the CEA executable for Windows only. On Linux and macOS
compile the Fortran source that comes with it once, and tell CEA_Wrap where
the executable is:

```sh
ASSETS=$(python -c "import CEA_Wrap, os; print(os.path.join(os.path.dirname(CEA_Wrap.__file__), 'assets'))")
mkdir -p ~/.local/bin && cd "$ASSETS"
gfortran -O2 -std=legacy -o ~/.local/bin/Dan_CEA2 dan_cea2.f
export CEA_EXE_LOCATION=~/.local/bin/Dan_CEA2    # add to ~/.bashrc
```

## Usage

```sh
ceatofoam design examples/rp1-lox-500psi.yaml          # print the design point
ceatofoam case   examples/rp1-lox-500psi.yaml cases/rp1 # write the case
cd cases/rp1 && ./Allrun                                # blockMesh ... rhoCentralFoam
ceatofoam post   cases/rp1                              # thrust, mass flow -> results.json
ceatofoam sweep  --fuel RP-1 --oxidizer "O2(L)" --pressure "500 psi" --plot of.png
ceatofoam interactive cases/mycase                     # console dialogue of version 1
```

`ceatofoam case` also writes `ceatofoam.json` into the case: the input,
the CEA results, the frozen mixture, the fit errors, the design point, the
contour, the mesh and the versions of every component (including the
source revision of ceaToFoam), so a case can always be traced back.

### Input

Quantities are SI numbers or strings with units (`"500 psi"`, `"2250 lbf"`,
`"40 in"`, `"15 deg"`).

```yaml
propellants: {fuel: RP-1, oxidizer: O2(L)}   # fuel_temperature, oxidizer_temperature optional
chamber:     {pressure: 500 psi, mixture_ratio: 2.4}
design:
  thrust: 500 lbf
  altitude: 0 ft          # or ambient_pressure and ambient_temperature
  sizing: frozen          # frozen | equilibrium | constant-gamma (gamma: gammas | frozen | <number>)
geometry:
  characteristic_length: 40 in
  convergent_half_angle: 30 deg
  divergent_half_angle: 15 deg
  contraction_ratio: auto         # 8.0 Dt[cm]^-0.6 + 1.25 (Humble et al. 1995)
  throat_upstream_radius: 1.5     # times the throat radius
  throat_downstream_radius: 1.5
  chamber_fillet_radius: 1.0
mesh:
  radial_cells: 40      # axis to wall
  cell_aspect: 1.0      # axial/radial cell size
  scale: 1.0            # multiplies every count (mesh studies)
  include_cylinder: false
  plume_length: 0       # exit diameters; 0: domain ends at the exit
  plume_radius: 3       # exit radii
solver:
  model: kOmegaSST      # inviscid | laminar | kOmegaSST
  iterations: 20000     # local time stepping: one time unit per iteration
  max_co: 0.2
```

A fixed-composition gas with a given chamber temperature replaces
`propellants` for cold-flow cases, and the geometry can be given
explicitly (`throat_radius`, `area_ratio`, `contraction_ratio`,
`chamber_length`); see the docstring of `ceatofoam/build.py`.

### The OpenFOAM case

* Thermophysics: `hePsiThermo`, `pureMixture`, `sutherland`, `janaf`,
  `perfectGas`, `sensibleInternalEnergy`. One pseudo-species: the frozen
  chamber mixture, which also fills the ambient region.
* Domain: from the end of the cylindrical chamber to the nozzle exit. In
  a steady frozen flow the cylinder is a duct at constant total
  conditions, and its low Mach number slows the density-based solver;
  `include_cylinder: true` meshes it anyway, and `plume_length` > 0 adds a
  plume region behind the exit (for flow pictures; thrust and mass flow
  need only the exit plane).
* Mesh: wedge of 5 degrees, symmetric about the x-y plane; patches `inlet`,
  `nozzle` (wall), `outlet` (the exit plane), `front`, `back`; with a plume
  also `base` (wall) and `farfield`.
* Boundary conditions: total pressure and temperature at the inlet,
  adiabatic walls, extrapolation at the supersonic exit; total-pressure
  entrainment at the plume farfield.
* Solution: `rhoCentralFoam` with local time stepping from the 1-D frozen
  solution. Function objects record the mass flow at the inlet, the throat
  and the exit plane, the exit-plane momentum and pressure integrals, the
  wall force, the extremes of p and T, and axis and wall profiles.
  `ceatofoam post` averages over the last quarter of the iterations and
  reports the thrust from the exit plane and, as a check of the momentum
  budget, from the wall force, with the mass balance and the iterative
  standard deviation. Judge convergence by the mass balance.

## Tests

```sh
python -m unittest discover -s tests
```

Tests that need the CEA executable or OpenFOAM are skipped when these are
not available.

## Changes from version 1.x

Version 2 is a rewrite of the three scripts of version 1 (tag `v1.0`).
It corrects the following:

* the nozzle was sized with CEA's equilibrium gamma_s and a molecular
  weight averaged over chamber, throat and exit, while the CFD solved the
  frozen chamber mixture with its own gamma; the nozzle was therefore
  over-expanded in the CFD and the exit pressure about 16 % below ambient;
* angles were converted to radians and then treated as degrees in the mesh,
  which collapsed the throat arcs and changed the convergent angle;
* the two wedge faces were rotated from already-rotated coordinates, which
  made the wedge asymmetric (checkMesh failed);
* the enthalpy and entropy constants of the janaf fit were evaluated at the
  chamber state instead of the reference state;
* the viscosity was that of air;
* the atmosphere-layer boundaries were in feet instead of metres;
* the mesh spacing and the plume domain did not scale with the engine.

## Limitations

The composition is frozen at the chamber (no finite-rate or shifting
chemistry in the nozzle), the gas is ideal, the wall is adiabatic, the
contour is conical, and the ambient gas is the exhaust mixture at ambient
temperature and pressure.

## License

GPL-3.0. Authors: Luke Cerda and Ehsan Madadi-Kandjani, California State
University Long Beach.
