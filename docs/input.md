# Input reference

The input is a YAML file (or, from Python, a dictionary) with the sections
below. Quantities are SI numbers or strings with units:

| Dimension | Units |
| --- | --- |
| pressure | `Pa`, `kPa`, `MPa`, `bar`, `atm`, `psi`, `psia` |
| force | `N`, `kN`, `MN`, `lbf` |
| length | `m`, `cm`, `mm`, `in`, `ft` |
| temperature | `K` |
| angle | `deg` |

A number without a unit is SI (angles in degrees). A unit of the wrong
dimension is an error.

## Complete example

```yaml
propellants:
  fuel: RP-1
  oxidizer: O2(L)
  # fuel_temperature: 298.15 K       # optional; CEA reference temperature otherwise
  # oxidizer_temperature: 90.17 K
chamber:
  pressure: 500 psi
  mixture_ratio: 2.4
design:
  thrust: 500 lbf
  altitude: 0 ft                     # or ambient_pressure and ambient_temperature
  sizing: frozen                     # frozen | equilibrium | constant-gamma
geometry:
  characteristic_length: 40 in
  convergent_half_angle: 30 deg
  divergent_half_angle: 15 deg
  contraction_ratio: auto
  throat_upstream_radius: 1.5        # times the throat radius
  throat_downstream_radius: 1.5
  chamber_fillet_radius: 1.0
mesh:
  radial_cells: 40
  cell_aspect: 1.0
  scale: 1.0
  wedge_angle: 5.0
  growth: 1.05
  include_cylinder: false
  plume_length: 0
  plume_radius: 3
solver:
  model: kOmegaSST
  iterations: 20000
  write_interval: 2000
  max_co: 0.2
  n_procs: 1
```

Every key except `propellants` (or `gas`), `chamber` and `design.thrust` has
the default shown.

## `propellants`

| Key | Default | Meaning |
| --- | --- | --- |
| `fuel` | required | CEA species name, e.g. `RP-1`, `CH4(L)`, `H2(L)`, `CH6N2(L)` |
| `oxidizer` | required | CEA species name, e.g. `O2(L)`, `N2O4(L)`, `O2` |
| `fuel_temperature` | CEA reference | inlet temperature of the fuel |
| `oxidizer_temperature` | CEA reference | inlet temperature of the oxidizer |

Liquid propellants enter at the reference temperature of their condensed
phase in the CEA library (90.17 K for O2(L), 111.643 K for CH4(L), 20.27 K for
H2(L), 298.15 K for RP-1), since CEA assigns their enthalpy at that
temperature.

### `gas` instead of `propellants`

For cold-flow cases (validation against air nozzles, for instance) a gas of
fixed composition replaces the propellants, and the chamber temperature is
given:

```yaml
gas:
  composition: {N2: 0.78084, O2: 0.20946, AR: 0.00934}   # mole fractions
  mechanism: air.yaml                                     # Cantera input file
chamber:
  pressure: 100 psi
  temperature: 294 K
```

## `chamber`

| Key | Default | Meaning |
| --- | --- | --- |
| `pressure` | required | chamber (stagnation) pressure \(p_c\) |
| `mixture_ratio` | required with `propellants` | oxidizer-to-fuel mass ratio O/F |
| `temperature` | required with `gas` | chamber temperature \(T_c\) |

## `design`

| Key | Default | Meaning |
| --- | --- | --- |
| `thrust` | required for a sized geometry | design thrust \(F_d\) at matched expansion |
| `altitude` | `0` | ambient state from the NASA Glenn Earth atmosphere model |
| `ambient_pressure`, `ambient_temperature` | — | ambient state given directly (temperature default 288.15 K) |
| `sizing` | `frozen` | sizing model, see below |
| `gamma` | `gammas` | constant-gamma model only: `gammas` (CEA equilibrium \(\gamma_s\) of the chamber), `frozen` (frozen \(\gamma\) at \(T_c\)) or a number |
| `molecular_weight` | `chamber` | constant-gamma model only: `chamber`, `station-average` (mean of the CEA chamber, throat and exit values) or a number [kg/kmol] |

The sizing models are described on the [Method](method.md#design-point-and-sizing-models)
page:

- `frozen` sizes the nozzle along the isentrope of the frozen chamber mixture
  with temperature-dependent \(c_p\). It is consistent with the simulation and
  is the recommended model.
- `equilibrium` takes \(I_{sp}\), \(c^*\) and the area ratio of CEA's
  shifting-equilibrium solution.
- `constant-gamma` uses the calorically perfect gas relations.

The last two are kept for comparison. In a frozen-flow simulation they give
nozzles that over-expand the gas (see [Verification and validation](validation.md#thermochemical-handoff)).

## `geometry`

The thrust chamber is sized from the design point unless `throat_radius` is
given.

| Key | Default | Meaning |
| --- | --- | --- |
| `characteristic_length` | `40 in` | \(L^*\); the chamber volume up to the throat is \(L^* A_t\) |
| `convergent_half_angle` | `30 deg` | \(\theta_c\) |
| `divergent_half_angle` | `15 deg` | \(\theta_d\) |
| `contraction_ratio` | `auto` | \(A_c/A_t\); `auto` uses \(8.0\,D_t^{-0.6} + 1.25\), \(D_t\) in cm (Humble et al., 1995) |
| `throat_upstream_radius` | `1.5` | radius of curvature upstream of the throat, in throat radii |
| `throat_downstream_radius` | `1.5` | radius of curvature downstream of the throat, in throat radii |
| `chamber_fillet_radius` | `1.0` | fillet between cylinder and convergent cone, in throat radii |

### Explicit geometry

To reproduce a given nozzle, give its dimensions instead of sizing it:

| Key | Meaning |
| --- | --- |
| `throat_radius` | \(r_t\) |
| `area_ratio` | \(A_e/A_t\) |
| `contraction_ratio` | \(A_c/A_t\) (a number) |
| `chamber_length` | optional length of the cylinder; otherwise \(L^*\) sets the volume |

The arc radii and cone angles above apply as well. The mass flow rate of the
one-dimensional frozen flow through this throat is then the reference for the
discharge coefficient.

## `mesh`

| Key | Default | Meaning |
| --- | --- | --- |
| `radial_cells` | `40` | cells between the axis and the wall, at every axial station |
| `cell_aspect` | `1.0` | axial over radial cell size |
| `scale` | `1.0` | multiplies every cell count (mesh studies) |
| `wedge_angle` | `5.0` | total angle of the wedge [deg] |
| `growth` | `1.05` | growth ratio of the cells in the plume region, from the exit-plane cell size |
| `include_cylinder` | `false` | mesh the cylindrical chamber as well |
| `plume_length` | `0` | length of a plume region behind the exit, in exit diameters; `0` ends the domain at the exit |
| `plume_radius` | `3` | radius of the plume region, in exit radii |

The cell counts scale with the throat radius, so the same settings give the
same resolution for a 500 lbf engine and for a 1.5 Mlbf engine.

## `solver`

| Key | Default | Meaning |
| --- | --- | --- |
| `model` | `kOmegaSST` | `inviscid` (slip walls), `laminar` or `kOmegaSST` (wall functions) |
| `iterations` | `20000` | iterations of local time stepping (one time unit each) |
| `write_interval` | `2000` | iterations between written fields |
| `max_co` | `0.2` | Courant number of the local time step |
| `limiter` | `vanLeer` | limiter of the reconstruction of \(\rho\), \(U\) and \(T\) |
| `lts_smoothing` | `0.02` | `rDeltaTSmoothingCoeff` |
| `turbulence_intensity` | `0.05` | at the inlet |
| `mixing_length` | `0.07` | at the inlet, in chamber diameters |
| `n_procs` | `1` | processors; above 1 `Allrun` decomposes and runs in parallel |
| `inlet_velocity_bc` | `zeroGradient` | or `pressureInletOutletVelocity` |

## Thermodynamic data

| Key | Default | Meaning |
| --- | --- | --- |
| `mechanism` (top level) | `gri30_highT.yaml` | Cantera data set in which the frozen chamber composition is evaluated |

Product species that the data set does not contain are dropped and the
remaining mole fractions renormalised. The dropped fraction is printed and
recorded in `ceatofoam.json`; for RP-1, CH₄ and H₂ with LOX at 500 and
1500 psia it does not exceed 20 ppm.
