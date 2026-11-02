# Wall heat flux

## Options

| Key | Meaning |
| --- | --- |
| `solver.wall_temperature` | wall temperature: a value [K], or a table of `[x - x_throat [m], T [K]]` rows (stock `fixedProfile`); without it the wall is adiabatic |
| `solver.wall_treatment` | `wall-function` or `resolved` |
| `solver.wall_output` | records the wall data on an adiabatic wall too |
| `mesh.first_cell_height` | height of the wall cell at the throat; it adds a graded layer |
| `mesh.layer_fraction` | graded outer fraction of the radius (default 0.25) |

The first cell scales with the local radius, so it is larger in the chamber
and at the exit than at the throat. Check \(y^+\) along the wall in the result.

```yaml
mesh: {radial_cells: 40, first_cell_height: 1.5e-6}
solver: {model: kOmegaSST, wall_treatment: resolved, wall_temperature: 500,
         lts_smoothing: 0.12, iterations: 60000, write_interval: 10000}
```

## Convergence

A wall layer needs two changes from the defaults. `lts_smoothing` must exceed the
growth of the layer: with the default 0.02 the thin wall cells hold back the time
step of the whole domain. And the wall heat flux settles far more slowly than the
mass flow: on the test nozzle the wall heat rate reached its plateau at about
50 000 iterations, and the coefficient at 20 000 was 26 % above its value at
40 000. Judge a run by `heat_rate_std` and `heat_rate_drift` in `results.json`.

On a resolved cooled wall the local flow does not become steady, so values along
the wall are averaged over the last quarter of the run, with their standard
deviation in `q_std`. In that case the inlet mass flow stayed 2 % above the throat
and exit flows; the cause is not established.

## Output

`ceatofoam post` writes `profile_wall_heat.csv` (heat flux into the wall,
\(y^+\), shear stress, temperature and pressure along the wall) and adds the
throat flux, the peak flux, the total heat rate and its convergence measures to
`results.json`.

`ceatofoam heat CASE_DIR` then writes `profile_bartz.csv`, which puts the
simulated flux next to the Bartz correlation evaluated with the frozen and,
for propellant cases, the equilibrium chamber gas.

## Correlation

Equations (4-10), (4-10-a), (4-13) and (4-14) of Huzel and Huang (1992),
pp. 85-86, in SI units, with \(r = Pr^{0.33}\) and the Mach number of the 1-D
expansion at the local area.

## Limits

The simulation is a frozen gas, so it does not contain the recombination near
a cold wall, which Huzel and Huang name as a reason for measured fluxes above
the calculated ones. OpenFOAM's Sutherland transport fixes the Prandtl number
through its Eucken form; compare it with the mixture value.
