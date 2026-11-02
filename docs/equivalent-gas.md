# Equivalent gas

With `design.gas: equilibrium` or `kinetic` the simulation runs a single
thermally perfect gas whose isentrope is the expansion of the reacting
mixture. `rhoCentralFoam`, `janaf` and `perfectGas` are used unchanged: there
is no species transport and no chemistry at run time.

## Idea

Along a streamline of an inviscid adiabatic flow,

\[
\mathrm{d}e = \frac{p}{\rho^2}\,\mathrm{d}\rho ,
\]

with \(e\) the internal energy including the heats of formation. This holds
whatever the chemistry does, so an expansion is a curve \(p(\rho)\),
\(e(\rho)\). With the gas constant \(R\) of the chamber define

\[
\tau = \frac{p}{\rho R}, \qquad e = e_{eq}(\tau).
\]

The gas \(p = \rho R \tau\), \(e = e_{eq}(\tau)\) has the same first law along
the same curve, so its isentrope through the chamber is the expansion, and its
speed of sound is \(\mathrm{d}p/\mathrm{d}\rho\) along it. At the chamber
\(\gamma\) of this gas equals CEA's \(\gamma_s\).

## Paths

- `equilibrium`: Cantera equilibrium at the chamber entropy.
- `kinetic`: Cantera rates through the contour of the case, with shifting
  equilibrium up to the inlet station and the pressure as marching variable,
  so that the sonic point is not singular. It depends on the size of the engine.

Both use the `gri30_highT.yaml` data. The rates of GRI-Mech 3.0 were not
developed for rocket chamber pressures; treat kinetic results accordingly.

## Fit

\(c_{p,eq} = c_{v,eq} + R\) is fitted by the two janaf quartics, both placed on
the \(\tau\) range of the path. `ceatofoam.json` records the fit error, the map
from \(\tau\) to the gas temperature and the 1-D expansion of the fitted gas,
which should reproduce the design.

## Limits

- **\(\tau\) is not the temperature.** The `T` field of the case is \(\tau\).
  `ceatofoam post` maps the exit value back and reports both.
- **Exact on the path only.** Strong shocks and wall layers leave the path.
- **Transport.** OpenFOAM's Sutherland transport derives the conductivity from
  \(c_v\), which here carries the heat of the reactions: viscous runs are
  estimates, and wall heat flux must not be taken from this gas.
- **One gas per chamber state**, as with the frozen mixture.
