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

"""OpenFOAM case for rhoCentralFoam (OpenFOAM.com, tested with v2606).

Thermophysics: hePsiThermo, pureMixture, sutherland transport, janaf
thermodynamics, perfectGas, sensibleInternalEnergy, with the coefficients
fitted to the frozen chamber mixture. The ambient gas is the same mixture at
ambient temperature and pressure (single-species model).

Flow models: ``inviscid`` (Sutherland As = 0, which rhoCentralFoam treats as
inviscid; slip walls), ``laminar`` and ``kOmegaSST`` (wall functions).

The case runs to steady state with local time stepping (localEuler; time is
the iteration count) from the 1-D frozen solution inside the chamber
(setFields), a column of exit-state gas behind the exit and ambient gas at
rest elsewhere. By default the domain ends at the nozzle exit (supersonic
outflow); with a plume region the exit plane is the faceZone ``exit``
(topoSet). Exit-plane integrals, the injector-face and throat mass flow and
the wall force on the patch ``nozzle`` are recorded by function objects;
ceatofoam.post turns them into thrust and mass flow.
"""

from dataclasses import dataclass
import math
import os
import stat

from .foamfile import fmt, write


@dataclass
class SolverSettings:
    model: str = "kOmegaSST"      # inviscid | laminar | kOmegaSST
    iterations: int = 20000
    write_interval: int = 2000
    max_co: float = 0.2
    limiter: str = "vanLeer"      # reconstruct() limiter for rho, U, T
    lts_smoothing: float = 0.02   # rDeltaTSmoothingCoeff (solver default)
    turbulence_intensity: float = 0.05
    mixing_length: float = 0.07   # times the chamber diameter
    n_procs: int = 1
    inlet_velocity_bc: str = "zeroGradient"   # or pressureInletOutletVelocity


@dataclass
class CaseConditions:
    """Thermodynamic and boundary data of one case (SI)."""

    pc: float
    tc: float
    p_amb: float
    t_amb: float
    mw: float
    gamma_c: float
    gamma_amb: float
    janaf: object
    sutherland: object
    x_throat: float
    x_exit: float
    rc: float
    re: float
    inlet_velocity: float
    initial_slices: list
    x_end: float
    wall_points: list
    initial_jet: tuple       # (p, T, u) of the 1-D exit state
    plume: bool              # mesh has the plume region behind the exit
    x_inlet: float           # axial position of the inlet patch


def _thermophysical(c, model):
    lo, hi = c.janaf.low, c.janaf.high
    As = 0.0 if model == "inviscid" else c.sutherland.As
    Ts = c.sutherland.Ts
    coeffs = lambda a: " ".join(fmt(v) for v in a)
    return f"""thermoType
{{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       sutherland;
    thermo          janaf;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;
}}

// Frozen chamber mixture. janaf coefficients and Sutherland constants are
// fitted by ceaToFoam (see ceatofoam.json in this case for the fit errors).
mixture
{{
    specie
    {{
        molWeight       {fmt(c.mw)};
    }}
    thermodynamics
    {{
        Tlow            {fmt(c.janaf.t_low)};
        Thigh           {fmt(c.janaf.t_high)};
        Tcommon         {fmt(c.janaf.t_common)};
        highCpCoeffs    ({coeffs(hi)});
        lowCpCoeffs     ({coeffs(lo)});
    }}
    transport
    {{
        As              {fmt(As)};
        Ts              {fmt(Ts)};
    }}
}}
"""


def _turbulence_properties(model):
    if model in ("inviscid", "laminar"):
        return "simulationType  laminar;\n"
    return f"""simulationType  RAS;

RAS
{{
    RASModel        {model};
    turbulence      on;
    printCoeffs     on;
}}
"""


_WEDGES = """    front
    {
        type            wedge;
    }
    back
    {
        type            wedge;
    }
"""


def _field(dims, internal, patches):
    return (
        f"dimensions      {dims};\n\n"
        f"internalField   uniform {internal};\n\n"
        "boundaryField\n{\n" + patches + _WEDGES + "}\n"
    )


def _p(c):
    return _field("[1 -1 -2 0 0 0 0]", fmt(c.p_amb), f"""    inlet
    {{
        type            totalPressure;
        p0              uniform {fmt(c.pc)};
        psi             thermo:psi;
        gamma           {fmt(c.gamma_c)};
        value           uniform {fmt(c.pc)};
    }}
    "(nozzle|base)"
    {{
        type            zeroGradient;
    }}
    // Supersonic outflow at the exit plane (no plume) or through the jet
    // core (plume). waveTransmissive failed at the plume outlet with local
    // time stepping in tests; extrapolation is stable.
    outlet
    {{
        type            zeroGradient;
    }}
    farfield
    {{
        type            totalPressure;
        p0              uniform {fmt(c.p_amb)};
        psi             thermo:psi;
        gamma           {fmt(c.gamma_amb)};
        value           uniform {fmt(c.p_amb)};
    }}
""")


def _T(c):
    # Initial face value at the outlet: the 1-D exit state when the domain
    # ends at the exit (setFields sets cells only).
    t_out = c.t_amb if c.plume else c.initial_jet[1]
    return _field("[0 0 0 1 0 0 0]", fmt(c.t_amb), f"""    inlet
    {{
        type            totalTemperature;
        T0              uniform {fmt(c.tc)};
        psi             thermo:psi;
        gamma           {fmt(c.gamma_c)};
        value           uniform {fmt(c.tc)};
    }}
    "(nozzle|base)"
    {{
        type            zeroGradient;
    }}
    outlet
    {{
        type            inletOutlet;
        inletValue      uniform {fmt(c.t_amb)};
        value           uniform {fmt(t_out)};
    }}
    farfield
    {{
        type            inletOutlet;
        inletValue      uniform {fmt(c.t_amb)};
        value           uniform {fmt(c.t_amb)};
    }}
""")


def _U(c, model, inlet_bc):
    wall = "slip" if model == "inviscid" else "noSlip"
    u_out = 0.0 if c.plume else c.initial_jet[2]
    return _field("[0 1 -1 0 0 0 0]", "(0 0 0)", f"""    inlet
    {{
        type            {inlet_bc};
        value           uniform (0 0 0);
    }}
    "(nozzle|base)"
    {{
        type            {wall};
    }}
    outlet
    {{
        type            inletOutlet;
        inletValue      uniform (0 0 0);
        value           uniform ({fmt(u_out)} 0 0);
    }}
    farfield
    {{
        type            pressureInletOutletVelocity;
        value           uniform (0 0 0);
    }}
""")


def _turbulence_fields(c, s):
    """k, omega, nut, alphat for kOmegaSST with wall functions."""
    k_in = 1.5*(s.turbulence_intensity*c.inlet_velocity)**2
    length = s.mixing_length*2*c.rc
    omega_in = math.sqrt(k_in)/(0.09**0.25*length)
    # Still ambient gas: 1 % of the inlet k with the inlet omega.
    k_amb, omega_amb = 1.0e-2*k_in, omega_in
    k = _field("[0 2 -2 0 0 0 0]", fmt(k_amb), f"""    inlet
    {{
        type            turbulentIntensityKineticEnergyInlet;
        intensity       {fmt(s.turbulence_intensity)};
        value           uniform {fmt(k_in)};
    }}
    "(nozzle|base)"
    {{
        type            kqRWallFunction;
        value           uniform {fmt(k_in)};
    }}
    outlet
    {{
        type            inletOutlet;
        inletValue      uniform {fmt(k_amb)};
        value           uniform {fmt(k_amb)};
    }}
    farfield
    {{
        type            inletOutlet;
        inletValue      uniform {fmt(k_amb)};
        value           uniform {fmt(k_amb)};
    }}
""")
    omega = _field("[0 0 -1 0 0 0 0]", fmt(omega_amb), f"""    inlet
    {{
        type            turbulentMixingLengthFrequencyInlet;
        mixingLength    {fmt(length)};
        value           uniform {fmt(omega_in)};
    }}
    "(nozzle|base)"
    {{
        type            omegaWallFunction;
        value           uniform {fmt(omega_in)};
    }}
    outlet
    {{
        type            inletOutlet;
        inletValue      uniform {fmt(omega_amb)};
        value           uniform {fmt(omega_amb)};
    }}
    farfield
    {{
        type            inletOutlet;
        inletValue      uniform {fmt(omega_amb)};
        value           uniform {fmt(omega_amb)};
    }}
""")
    calc = """    inlet
    {
        type            calculated;
        value           uniform 0;
    }
    outlet
    {
        type            calculated;
        value           uniform 0;
    }
    farfield
    {
        type            calculated;
        value           uniform 0;
    }
"""
    nut = _field("[0 2 -1 0 0 0 0]", "0", calc + """    "(nozzle|base)"
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
""")
    alphat = _field("[1 -1 -1 0 0 0 0]", "0", calc + """    "(nozzle|base)"
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
""")
    return {"k": k, "omega": omega, "nut": nut, "alphat": alphat}


def _control_dict(c, s):
    wall = "\n".join(
        f"                    ({fmt(x)} {fmt(r)} 0)" for x, r in c.wall_points
    )
    # The exit plane is the faceZone exit (topoSet) inside a plume domain,
    # or the outlet patch when the domain ends at the exit.
    exit_region, exit_name = (
        ("faceZone", "exit") if c.plume else ("patch", "outlet")
    )
    return f"""application     rhoCentralFoam;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {s.iterations};

// Local time stepping (ddtSchemes localEuler): one time unit per iteration.
deltaT          1;
maxCo           {fmt(s.max_co)};
maxDeltaT       1;
rDeltaTSmoothingCoeff {fmt(s.lts_smoothing)};

writeControl    timeStep;
writeInterval   {s.write_interval};
purgeWrite      2;
writeFormat     ascii;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable true;

functions
{{
    // Extremes of T and p with their locations (diagnostics).
    minMax
    {{
        type            fieldMinMax;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   10;
        log             false;
        location        true;
        mode            magnitude;
        fields          (T p);
    }}

    MachNo
    {{
        type            MachNo;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
    }}

    // Mass flow through the injector face (negative: inflow).
    inletFlux
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   10;
        log             false;
        writeFields     false;
        regionType      patch;
        name            inlet;
        operation       sum;
        fields          (phi);
        writeArea       true;
    }}
    inletPressure
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   10;
        log             false;
        writeFields     false;
        regionType      patch;
        name            inlet;
        operation       areaIntegrate;
        fields          (p);
    }}
    inletMomentum
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   10;
        log             false;
        writeFields     false;
        regionType      patch;
        name            inlet;
        operation       weightedSum;
        weightField     phi;
        fields          (U);
    }}

    // Nozzle exit plane (faceZone exit, oriented +x): mass flow, momentum
    // flux sum(phi U), pressure integral and area.
    exitFlux
    {{
        $inletFlux;
        regionType      {exit_region};
        name            {exit_name};
    }}
    exitPressure
    {{
        $inletPressure;
        regionType      {exit_region};
        name            {exit_name};
    }}
    exitMomentum
    {{
        $inletMomentum;
        regionType      {exit_region};
        name            {exit_name};
    }}
    exitAverage
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   10;
        log             false;
        writeFields     false;
        regionType      {exit_region};
        name            {exit_name};
        operation       areaAverage;
        fields          (p T Ma);
    }}

    // Profiles at every write: along the axis (just off it, inside the
    // first cell row) and on the chamber wall (nearest wall faces).
    profiles
    {{
        type            sets;
        libs            (sampling);
        writeControl    writeTime;
        setFormat       raw;
        interpolationScheme cell;
        fields          (p T Ma U);
        sets
        {{
            axis
            {{
                type        uniform;
                axis        x;
                start       ({fmt(c.x_inlet)} {fmt(1e-3*c.re)} 0);
                end         ({fmt(c.x_end)} {fmt(1e-3*c.re)} 0);
                nPoints     400;
            }}
            wall
            {{
                type        patchCloud;
                axis        x;
                patches     (nozzle);
                maxDistance {fmt(0.5*c.re)};
                points
                (
{wall}
                );
            }}
        }}
    }}

    // Mass flow through the throat plane (cutting plane, integral of
    // rhoU . dA): with inlet and exit it shows where the flow is unsteady.
    throatFlux
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   10;
        log             false;
        writeFields     false;
        regionType      sampledSurface;
        name            throat;
        sampledSurfaceDict
        {{
            type            cuttingPlane;
            planeType       pointAndNormal;
            pointAndNormalDict
            {{
                point       ({fmt(c.x_throat)} 0 0);
                normal      (1 0 0);
            }}
            interpolate     false;
        }}
        operation       areaNormalIntegrate;
        fields          (rhoU);
        writeArea       true;
    }}

    // Pressure (relative to ambient) and viscous force of the gas on the
    // chamber wall.
    wallForce
    {{
        type            forces;
        libs            (forces);
        writeControl    timeStep;
        writeInterval   10;
        log             false;
        patches         (nozzle);
        rho             rho;
        pRef            {fmt(c.p_amb)};
        CofR            (0 0 0);
    }}
}}
"""


def _fv_schemes(s):
    lim = s.limiter
    # Limited schemes have a vector form (vanLeerV, MinmodV); upwind has not.
    lim_v = lim if lim in ("upwind", "linear") else f"{lim}V"
    turb = ""
    if s.model == "kOmegaSST":
        turb = """    div(phi,k)      Gauss upwind;
    div(phi,omega)  Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
"""
    return f"""fluxScheme      Kurganov;

ddtSchemes
{{
    default         localEuler;
}}

gradSchemes
{{
    default         Gauss linear;
}}

divSchemes
{{
    default         none;
    div(tauMC)      Gauss linear;
{turb}}}

laplacianSchemes
{{
    default         Gauss linear corrected;
}}

interpolationSchemes
{{
    default         linear;
    reconstruct(rho) {lim};
    reconstruct(U)  {lim_v};
    reconstruct(T)  {lim};
}}

snGradSchemes
{{
    default         corrected;
}}

wallDist
{{
    method          meshWave;
}}
"""


_FV_SOLUTION = """solvers
{
    "(rho|rhoU|rhoE)"
    {
        solver          diagonal;
    }

    U
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        nSweeps         2;
        tolerance       1e-10;
        relTol          0;
    }

    e
    {
        $U;
        tolerance       1e-10;
        relTol          0;
    }

    "(k|omega)"
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        nSweeps         2;
        tolerance       1e-10;
        relTol          0;
    }
}
"""


def _top_set(c):
    return f"""actions
(
    {{
        name            exit;
        type            faceZoneSet;
        action          new;
        source          planeToFaceZone;
        point           ({fmt(c.x_exit)} 0 0);
        normal          (1 0 0);
        option          all;
    }}
);
"""


def _set_fields(c):
    """Initial field: the 1-D frozen solution in axial slices inside the
    chamber (``c.initial_slices``: list of (x0, x1, p, T, u)), ambient gas
    at rest in the plume."""
    big = 10*max(c.rc, c.re)
    regions = []
    for x0, x1, p, T, u in c.initial_slices:
        regions.append(f"""    boxToCell
    {{
        box ({fmt(x0)} {fmt(-big)} {fmt(-big)}) ({fmt(x1)} {fmt(big)} {fmt(big)});
        fieldValues
        (
            volScalarFieldValue p {fmt(p)}
            volScalarFieldValue T {fmt(T)}
            volVectorFieldValue U ({fmt(u)} 0 0)
        );
    }}""")
    # The jet is started as a column of exit-state gas downstream of the
    # exit: an impulsive start into still gas rolls up a starting vortex
    # whose low-pressure core reached the farfield boundary and failed
    # there (tests with a 3 r_e farfield, iteration ~3900).
    pj, tj, uj = c.initial_jet
    if c.plume:
        regions.append(f"""    cylinderToCell
    {{
        point1 ({fmt(c.x_exit)} 0 0);
        point2 ({fmt(c.x_end + 1.0)} 0 0);
        radius {fmt(c.re)};
        fieldValues
        (
            volScalarFieldValue p {fmt(pj)}
            volScalarFieldValue T {fmt(tj)}
            volVectorFieldValue U ({fmt(uj)} 0 0)
        );
    }}""")
    return f"""defaultFieldValues
(
    volScalarFieldValue p {fmt(c.p_amb)}
    volScalarFieldValue T {fmt(c.t_amb)}
    volVectorFieldValue U (0 0 0)
);

regions
(
{chr(10).join(regions)}
);
"""


def _decompose(n):
    return f"numberOfSubdomains {n};\n\nmethod          scotch;\n"


_ALLRUN = """#!/bin/sh
cd "${0%/*}" || exit 1
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions

restore0Dir
runApplication blockMesh
runApplication checkMesh
runApplication topoSet
runApplication setFields
if [ "$(getNumberOfProcessors)" -gt 1 ]
then
    runApplication decomposePar
    runParallel $(getApplication)
    runApplication reconstructPar -latestTime
else
    runApplication $(getApplication)
fi
"""

_ALLCLEAN = """#!/bin/sh
cd "${0%/*}" || exit 1
. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions
cleanCase0
"""


def write_case(case, conditions, solver):
    """Write every dictionary except system/blockMeshDict."""
    c, s = conditions, solver
    write(case, "constant/thermophysicalProperties", "dictionary",
          _thermophysical(c, s.model))
    write(case, "constant/turbulenceProperties", "dictionary",
          _turbulence_properties(s.model))
    write(case, "0.orig/p", "volScalarField", _p(c))
    write(case, "0.orig/T", "volScalarField", _T(c))
    write(case, "0.orig/U", "volVectorField",
          _U(c, s.model, s.inlet_velocity_bc))
    if s.model == "kOmegaSST":
        for name, body in _turbulence_fields(c, s).items():
            write(case, f"0.orig/{name}", "volScalarField", body)
    write(case, "system/controlDict", "dictionary", _control_dict(c, s))
    write(case, "system/fvSchemes", "dictionary", _fv_schemes(s))
    write(case, "system/fvSolution", "dictionary", _FV_SOLUTION)
    write(case, "system/topoSetDict", "dictionary", _top_set(c))
    write(case, "system/setFieldsDict", "dictionary", _set_fields(c))
    write(case, "system/decomposeParDict", "dictionary",
          _decompose(s.n_procs))
    for name, text in (("Allrun", _ALLRUN), ("Allclean", _ALLCLEAN)):
        path = os.path.join(case, name)
        with open(path, "w") as f:
            f.write(text)
        os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP)
