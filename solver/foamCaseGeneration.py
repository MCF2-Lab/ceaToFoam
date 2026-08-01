
from ceaToFoam import a_sound, gamma, T1, Tamb, p0, P1, meanMolarMass, Tlow, Thigh, Tcommon, Hcof, Lcof, high_enthalpy_offset, high_entropy_offset, low_enthalpy_offset, low_entropy_offset, L_cylindrical, Rc, Rt, R_throat, Re, Lconv, Ldiv, convergent_half_angle, divergent_half_angle
import os
import shutil
import numpy as np
import math
import matplotlib.pyplot as plt
import io

##Asssuming a Mach number of 0.1 in the chamber for the velocity input and finding the velocity at the inlet by multiplying the speed of sound by this Mach.
VELOCITY_INPUT = 0.1*a_sound
U = "U.txt"
U_text = r'''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{
    inlet
    {
        type            zeroGradient;
    }

    outlet
    {
        type            waveTransmissive;
        gamma           GAMMA_VAL;
    }

    nozzle
    {
        type            noSlip;
    }
    asym2
    {
        type            wedge;
        
    }
    asym1
    {
        type            wedge;
        
    }
    outlet_r
    {
        type            waveTransmissive;
        gamma           GAMMA_VAL;
    }
}


// ************************************************************************* //'''

GAMMA_INPUT = str(gamma)
U_text = U_text.replace("XVEL", str(VELOCITY_INPUT)) ##replaces the XVEL in the U file with the user inputted value for the velocity in the x-direction
U_text = U_text.replace("GAMMA_VAL", GAMMA_INPUT) ##replaces the GAMMA_VAL in the U file with the user inputted value for the gamma

U_file = open(U, 'w')
U_file.write(U_text)
U_file.close()


TEMPERATURE_INPUT = str(T1)
TEMPERATURE_AMBIENT = str(Tamb + 273.15) ##converts the ambient temperature from Celsius to Kelvin
T = "T.txt"
T_text = r'''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;
    location    "0";
    object      T;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   uniform T_AMBIENT;
boundaryField
{
    outlet
    {
        type            waveTransmissive;
        gamma           GAMMA_VAL;
        fieldInf        T_AMBIENT;
        lInf            10;
        value           uniform T_AMBIENT;
    }
    inlet
    {
        type         totalTemperature;
	    gamma        GAMMA_VAL;
        T0           uniform T_VAL;
	    value	     uniform T_VAL;
    }
    nozzle
    {
        type            zeroGradient;
    }
    asym2
    {
        type            wedge;
        
    }
    asym1
    {
        type            wedge;
        
    }

    outlet_r
    {
        type            waveTransmissive;
        gamma           GAMMA_VAL;
        fieldInf        T_AMBIENT;
        lInf            10;
        value           uniform T_AMBIENT;
    }
}


// ************************************************************************* //'''

T_text = T_text.replace("T_VAL", TEMPERATURE_INPUT) ##replaces the TVAL in the T file with the user inputted value for the temperature
T_text = T_text.replace("T_AMBIENT", TEMPERATURE_AMBIENT) ##replaces the T_AMBIENT in the T file with the user inputted value for the ambient temperature
T_text = T_text.replace("GAMMA_VAL", GAMMA_INPUT) ##replaces the GAMMA_VAL in the T file with the user inputted value for the gamma

T_file = open(T, 'w')
T_file.write(T_text)
T_file.close()





alphat = "alphat.txt"
alphat_text = r'''
/*--------------------------------*- C++ -*----------------------------------*| =========                 |                                                 |
| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \    /   O peration     | Version:  v2412                                 |
|   \  /    A nd           | Website:  www.openfoam.com                      |
|    \/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      alphat;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            calculated;
        value           uniform 0;
    }

    outlet
    {
        type            calculated;
        value           uniform 0;
    }

    nozzle
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }

    asym2
    {
        type            wedge;
        
    }
    asym1
    {
        type            wedge;
        
    }
    outlet_r
    {
        type            calculated;
        value           uniform 0;
    }

}


// ************************************************************************* //'''

alphat_file = open(alphat, 'w')
alphat_file.write(alphat_text)
alphat_file.close()


epsilon = "epsilon.txt"
epsilon_text = r'''
/*--------------------------------*- C++ -*----------------------------------*| =========                 |                                                 |
| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \    /   O peration     | Version:  v2412                                 |
|   \  /    A nd           | Website:  www.openfoam.com                      |
|    \/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      epsilon;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -3 0 0 0 0];

internalField   uniform 266000;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 266000;
    }

    outlet
    {
        type            inletOutlet;
        inletValue      uniform 266000;
        value           uniform 266000;
    }
    outlet_r
    {
        type            inletOutlet;
        inletValue      uniform 266000;
        value           uniform 266000;
    }
    nozzle
    {
        type            epsilonWallFunction;
        value           uniform 266000;
    }

    asym2
    {
        type            wedge;
        
    }
    asym1
    {
        type            wedge;
        
    }
}


// ************************************************************************* //'''

epsilon_file = open(epsilon, 'w')
epsilon_file.write(epsilon_text)
epsilon_file.close()



k = "k.txt"

k_text = r'''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 1000;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 1000;
    }

    outlet
    {
        type            inletOutlet;
	inletValue 			uniform 1000;
        value           uniform 1000;
    }

    nozzle
    {
        type            kqRWallFunction;
        value           uniform 1000;
    }

    asym2
    {
        type            wedge;
        
    }
    asym1
    {
        type            wedge;
        
    }
    outlet_r
    {
        type            inletOutlet;
        inletValue      uniform 1000;
        value           uniform 1000;
    }

}


// ************************************************************************* //'''

k_file = open(k, 'w')
k_file.write(k_text)
k_file.close()



nut = "nut.txt"

nut_text = r'''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            calculated;
        value           uniform 0;
    }

    outlet
    {
        type            calculated;
        value           uniform 0;
    }
    outlet_r
    {
        type            calculated;
        value           uniform 0;
    }

    nozzle
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    asym2
    {
        type            wedge;
        
    }
    asym1
    {
        type            wedge;
        
    }
}


// ************************************************************************* //'''

nut_file = open(nut, 'w')
nut_file.write(nut_text)
nut_file.close()


pressure = "p.txt"
ambient_pressure_str = str(p0)
 
pressure_text = r'''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;
    location    "0";
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform AMB_PRESSURE;


boundaryField
{
    outlet
    {
        type            waveTransmissive;
        field           p;
        gamma           GAMMA_VAL;
        fieldInf        AMB_PRESSURE;
        psi             thermo:psi;
        lInf            1;
        value           uniform AMB_PRESSURE;
    }
    inlet
    {
        type            totalPressure;
        gamma           GAMMA_VAL;
        psi             thermo:psi;
        p0              uniform P_VAL;
        value           uniform P_VAL;
    }
    nozzle
    {
        type            zeroGradient;
    }
    outlet_r
    {
        type            waveTransmissive;
	    field           p;
        gamma           GAMMA_VAL;
        fieldInf        AMB_PRESSURE;
	    psi             thermo:psi;
        lInf            1;
        value           uniform AMB_PRESSURE;
    }
    asym2
    {
        type            wedge;
        
    }
    asym1
    {
        type            wedge;
        
    }
}


// ************************************************************************* //'''

replacements = {"AMB_PRESSURE": ambient_pressure_str,
                "P_VAL": str(P1),
                "GAMMA_VAL": GAMMA_INPUT
}
for old, new in replacements.items(): ##this replaces the ambient pressure and chamber pressure specified by the user into the pressure file for OpenFOAM
     pressure_text = pressure_text.replace(old, new)
 
pressure_file = open(pressure, 'w')
pressure_file.write(pressure_text)
pressure_file.close()
 
thermophysicalProperties = "thermophysicalProperties.txt"

tpp_text = r'''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

thermoType
{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       sutherland;
    thermo          janaf;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;
}

mixture
{
    specie
    {
        nMoles          1;
        molWeight       MOLECULAR_WEIGHT;
    }
    thermodynamics
    {
        Tlow            T_LOW;
        Thigh           T_HIGH;
        Tcommon         T_COMMON;
        highCpCoeffs    ( Hcof0 Hcof1 Hcof2 Hcof3 Hcof4  high_enthalpy_offset high_entropy_offset );
        lowCpCoeffs     ( Lcof0 Lcof1 Lcof2 Lcof3 Lcof4  low_enthalpy_offset low_entropy_offset );
    }
    transport
    {
        As              1.67212e-06;
        Ts              170.672;
    }
}


// ************************************************************************* //'''

replacements_tpp = {"MOLECULAR_WEIGHT": str(meanMolarMass),
                    "T_LOW": str(Tlow),
                    "T_HIGH": str(Thigh),
                    "T_COMMON": str(Tcommon),
                    "Hcof0": str(Hcof[0]),
                    "Hcof1": str(Hcof[1]),
                    "Hcof2": str(Hcof[2]),
                    "Hcof3": str(Hcof[3]),
                    "Hcof4": str(Hcof[4]),
                    "high_enthalpy_offset": str(high_enthalpy_offset),
                    "high_entropy_offset": str(high_entropy_offset),
                    "Lcof0": str(Lcof[0]),
                    "Lcof1": str(Lcof[1]),
                    "Lcof2": str(Lcof[2]),
                    "Lcof3": str(Lcof[3]),
                    "Lcof4": str(Lcof[4]),
                    "low_enthalpy_offset": str(low_enthalpy_offset),
                    "low_entropy_offset": str(low_entropy_offset),
                    "P_VAL": str(P1)
}

print(replacements_tpp)



for old, new in replacements_tpp.items(): ##this replaces the ambient pressure and chamber pressure specified by the user into the pressure file for OpenFOAM
    tpp_text = tpp_text.replace(old, new)

thermo_file = open(thermophysicalProperties, 'w')
thermo_file.write(tpp_text)
thermo_file.close()


turbulenceProperties = "turbulenceProperties.txt"

turbProp_text = r'''
/*--------------------------------*- C++ -*----------------------------------*| =========                 |                                                 |
| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \    /   O peration     | Version:  v2412                                 |
|   \  /    A nd           | Website:  www.openfoam.com                      |
|    \/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
simulationType  RAS;

RAS
{
    RASModel        kEpsilon;

    turbulence      on;

    printCoeffs     on;
}

// ************************************************************************* //'''

turbulence_file = open(turbulenceProperties, 'w')
turbulence_file.write(turbProp_text)
turbulence_file.close()


control_dict = "controlDict.txt"

control_dict_text = r'''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     rhoCentralFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

functions
{
    MachNo
    {
        type            MachNo;
        libs            ("libfieldFunctionObjects.so");
        writeControl    runTime;
        writeInterval   5e-05;     
    }
    wallHeatFlux
    {
        type            wallHeatFlux;
        libs            ("libfieldFunctionObjects.so");
        writeControl    runTime;
        writeInterval   5e-05;
        patches        (nozzle);
    }
}

endTime         5.0e-02;

deltaT          2e-07;

writeControl    runTime;

writeInterval   5e-05;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

adjustTimeStep  yes;

maxCo           0.5;

//maxDeltaT       1;


// ************************************************************************* //'''


control_dict_file = open(control_dict, 'w')
control_dict_file.write(control_dict_text)
control_dict_file.close()



fv_schemes = "fvSchemes.txt"

fv_schemes_text = r'''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

fluxScheme      Kurganov;

ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(tauMC)      Gauss linear;
    div(phi,U)      Gauss limitedLinearV 1;
    div(phi,e)      Gauss limitedLinear 1;
    div(phid,p)     Gauss limitedLinear 1;
    div(phi,K)      Gauss limitedLinear 1;
    div(phiv,p)     Gauss limitedLinear 1;
    div(phi,k)      Gauss upwind;
    div(phi,epsilon) Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
    reconstruct(rho) vanAlbada;
    reconstruct(U)  vanAlbadaV;
    reconstruct(T)  vanAlbada;
}

snGradSchemes
{
    default         corrected;
}


// ************************************************************************* //'''


fv_schemes_file = open(fv_schemes, 'w')
fv_schemes_file.write(fv_schemes_text)
fv_schemes_file.close()



fv_solutions = "fvSolution.txt"

fv_solutions_text = r'''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
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
        relTol          0.0;
    }

    e
    {
        $U;
        tolerance       1e-10;
        relTol          0.0;
    }
    relaxationFactors
    {
    equations
    {
          rho      0.15;
          rhoU     0.15;
      rhoE     0.15;
    }
    }
    "(k|epsilon).*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-08;
        relTol          0;
    }
}

// ************************************************************************* //'''


fv_solutions_file = open(fv_solutions, 'w')
fv_solutions_file.write(fv_solutions_text)
fv_solutions_file.close()




decomposePar = "decomposeParDict.txt"

decomposePar_text = r'''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains 10;

method          scotch;

// ************************************************************************* //'''


decomposePar_file = open(decomposePar, 'w')
decomposePar_file.write(decomposePar_text)
decomposePar_file.close()

comb_r=Rc ##radius of combustion chamber
comb_l=L_cylindrical ##straight length of combustion chamber
comb_r1=0.001## curvature radius of combustion chamber tapering
comb_half_angle= convergent_half_angle##half angle of the combustion chamber tapering
throat_r=Rt ## radius at the nozzle
comb_full_l=L_cylindrical+Lconv ##length of combustion chamber till throat
throat_curv_r=R_throat ## curvature radius at the nozzle throat
div_half_angle_1=divergent_half_angle ##half angle of the diverging section
straight_1=0.001 ## first transitional straight section
div_r = 0  ##no bell curvature for a purely conical nozzle
straight_2=0.0014 ##straight section at the nozzle exit
div_half_angle_2=div_half_angle_1 ##half angle of the diverging section at the termination
total_l=(L_cylindrical+Lconv+Ldiv) ##total length of the engine
bell_r=Re ##radius of the nozzle bell at the rim
domain_r=10*bell_r ##radius of the area behind the nozzle for the plume
domain_ext=1.5*total_l+6 ##extra length for the plume
wedge_half_angle=2.5 ##half angle to create a wedge

points = np.zeros((31, 3)) ##main blockmesh vertices
aux_points = np.zeros((6, 3)) ##auxiliary points to specify arc segments

##points 0-8 are placed on the wedge axis

points[0]=(0,0,0)
points[1]=(comb_l,0,0)
points[2]=(comb_l+comb_r1*math.sin(math.pi*comb_half_angle/180),0,0)
points[3]=(comb_full_l-throat_curv_r*math.sin(math.pi*comb_half_angle/180),0,0)
points[4]=(comb_full_l+throat_curv_r*math.sin(math.pi*div_half_angle_1/180),0,0)
points[5]=(points[4,0]+straight_1*math.cos(math.pi*div_half_angle_1/180),0,0)
points[6]=(total_l-straight_2*math.cos(math.pi*div_half_angle_2/180),0,0)
points[7]=(total_l,0,0)
points[8]=(total_l+domain_ext,0,0)

for i in range(9,17):
    points[i,0]=points[i-9,0]
points[9,1]=comb_r
points[10,1]=comb_r
points[11,1]=comb_r-comb_r1*(1-math.cos(math.pi*comb_half_angle/180))
points[12,1]=throat_r+throat_curv_r*(1-math.cos(math.pi*comb_half_angle/180))
points[13,1]=throat_r+throat_curv_r*(1-math.cos(math.pi*div_half_angle_1/180))
points[14,1]=points[13,1]+straight_1*math.sin(math.pi*div_half_angle_1/180)
points[15,1]=bell_r-straight_2*math.sin(math.pi*div_half_angle_2/180)
points[16,1]=bell_r

points[17]=(total_l,domain_r,0)
points[18]=(total_l+domain_ext,domain_r,0)
points[19]=(total_l+domain_ext,bell_r,0)

#plt.style.use('fivethirtyeight')

#plt.plot(points[9:20,0], points[9:20,1])
#plt.plot(points[9:20,0], -1*points[9:20,1])
plt.plot([points[0,0],points[9,0]], [points[0,1],points[9,1]], 'b')
plt.plot(points[9:18,0], points[9:18,1], 'r')
plt.plot(points[17:19,0], points[17:19,1], 'k')
plt.plot([points[18,0],points[19,0],points[8,0]], [points[18,1],points[19,1],points[8,1]], 'm')


#plt.plot(points[0:9,0], points[0:9,1])

plt.title(f'geometry with colors representing boundary condition zones')
plt.xlim(left=-0.1)                ##set up lower y-axis limit at zero
plt.xlim(right=10)
                   ##set upper limit of the graph at 110% of maximum beam radius
plt.xlabel('Distance, [m]')
plt.ylabel('Radius, [m]')


#plt.savefig('Rothe_nozzle_profile.png', dpi=300)

##arc central point in combustion chamber
aux_points[0]=(comb_l+comb_r1*(math.sin(math.pi*comb_half_angle/360)),comb_r-comb_r1*(1-math.cos(math.pi*comb_half_angle/360)), 0)
##arc central point in the throat
aux_points[1]=(comb_full_l,throat_r, 0)

##arc central point in the diverging bell - only if a bell arc is requested
if div_r is not None and div_r > 0:
    midangle = (div_half_angle_1 + div_half_angle_2) / 2.0
    x_center = points[14, 0] + div_r * math.sin(math.pi * div_half_angle_1 / 180.0)
    y_center = points[14, 1] - div_r * math.cos(math.pi * div_half_angle_1 / 180.0)
    aux_points[2] = (x_center - div_r * math.sin(math.pi * midangle / 180.0),
                     y_center + div_r * math.cos(math.pi * midangle / 180.0), 0)
else:
    ##no bell arc for pure cone
    aux_points[2] = (math.nan, math.nan, 0)

#plt.scatter(aux_points[0:3,0], aux_points[0:3,1])
plt.show()

captured_output = io.StringIO()  ##Create a StringIO object to capture the output


## convert points to wedge

def rotate_points(points,angle):
    for point in points:
        point[1]=point[1]*math.cos(math.pi*angle/180)
        point[2]=point[1]*math.sin(math.pi*angle/180)
    return points

##generate main and aux wedge points by rotation
points[20:31]=rotate_points(points[9:20],-1*wedge_half_angle )
points[9:20]=rotate_points(points[9:20],wedge_half_angle )
aux_points[3:6]=rotate_points(aux_points[0:3],-1*wedge_half_angle )
aux_points[0:3]=rotate_points(aux_points[0:3],wedge_half_angle )

print("Blockmeshdict segment with produced vertices")
print("vertices")
print("(")
vertices_lines = []
for i, Point in enumerate(points):
    line = f"\t({Point[0]} {Point[1]} {Point[2]}) // {i} "
    print(line)                     ##keep console output
    vertices_lines.append(line)

print(");")

##now create the vertices_text used in replacements
vertices_text = "\n".join(vertices_lines)


axis_point_numbers=np.arange(0,9) ##poit numbers on the axis
wedge_point_numbers_1=np.arange(9,20) ## point numbers on the first wedge
wedge_point_numbers_2=np.arange(20,31) ##point numbers on the second wedge
##number of generated blocks
n_blocks=7

## define function that creates pseudo 4-point cross-sections for the main engine part
def contour(i):
    return [axis_point_numbers[i], wedge_point_numbers_2[i], wedge_point_numbers_1[i], axis_point_numbers[i]]

contours=[contour(i) for i in range(n_blocks+1)] ##define list of contours
blocks=[contours[i]+contours[i+1] for i in range(7)] ##define list of blocks

y_cells=[20 for i in range(n_blocks)] ##define list of cell numbers in y-direction  (uniform)
x_cells=[int((points[i+1,0]-points[i,0])/0.005+1) for i in range(n_blocks)] ##define list of cell numbers in x-direction  (every 5 mm)

##generate block section of blockmeshdict
print("blocks")
print("(")

blocks_text = []
for i,block in enumerate(blocks):
    
    block = (f'\t hex ({block[0]} {block[1]} {block[2]} {block[3]} {block[4]} {block[5]} {block[6]} {block[7]} ) ({y_cells[i]} 1 {x_cells[i]}) simpleGrading (1 1 1)')
    print(block)
    blocks_text.append(block)

    blocks_text2 = "\n".join(blocks_text)
        


##number of generated blocks

##generate asym1 surfaces for wedge BC
asym1=[[contours[i][0],contours[i][2],contours[i+1][2],contours[i+1][0]] for i in range(n_blocks) ]
print (f'asym1 patches')

asym_text22 = []
for patch in asym1:
    asym1_text = (f'\t ({patch[0]} {patch[1]} {patch[2]} {patch[3]})')
    print(asym1_text)
    asym_text22.append(asym1_text)
    asym1_text2 = "\n".join(asym_text22)

##generate asym2 surfaces for wedge BC

asym2=[[contours[i+1][0],contours[i+1][1],contours[i][1],contours[i][0]] for i in range(n_blocks) ]
print (f'asym2 patches')
asym_text = []
for patch in asym2:
    asym_text2 = (f'\t({patch[0]} {patch[1]} {patch[2]} {patch[3]})')
    print(asym_text2)
    asym_text.append(asym_text2)
    asym2_text = "\n".join(asym_text)
    
##generate nozzle surfaces for nozzle wall BC
##nozzle[i]=np.array((wedges[i,2],wedges[i,1], wedges[i+1,1], wedges[i+1,2]))
nozzle_text = []
nozzle=[[contours[i][2],contours[i][1],contours[i+1][1],contours[i+1][2]] for i in range(n_blocks) ]
print (f'nozzle patches')
for patch in nozzle:
    nozzle2 = (f'\t({patch[0]} {patch[1]} {patch[2]} {patch[3]})')
    print(nozzle2)
    nozzle_text.append(nozzle2)
    nozzle_text2 = "\n".join(nozzle_text)

##generate inlet surfaces for inlet BC
print (f'inlet patch')
inlet=contours[0]
inlet_text = []
inlet2 = (f'\t({inlet[0]} {inlet[1]} {inlet[2]} {inlet[3]} )')
print(inlet2)
inlet_text.append(inlet2)
inlet_text2 = "\n".join(inlet_text)


print(f'outlet_r patch')
print(
    """
           //(7 27 16 7 ) 
           // (8 30 19 8)//
            //(30 29 18 19)//
            (17 28 29 18) //
"""
)

##generate outlet surfaces for outlet BC

print (f'outlet patch')
outlet=contours[n_blocks]
print(f'({outlet[0]} {outlet[1]} {outlet[2]} {outlet[3]} )')
outlet_text = captured_output.getvalue()


print(f'arcs for the nozzle')

print("edges")
print("(")

edge1_text = []
edge2_text = []
edge3_text = []
edge4_text = []
edge5_text = []


edge1_2 = (f'\t arc 10 11 ({aux_points[0,0]} {aux_points[0,1]} {aux_points[0,2]})')
edge1_text.append(edge1_2)
edge12_text = "\n".join(edge1_text)


edge2_2 =(f'\t arc 12 13 ({aux_points[1,0]} {aux_points[1,1]} {aux_points[1,2]})')
edge2_text.append(edge2_2)
edge2_text2 = "\n".join(edge2_text)




if div_r is not None and div_r > 0:
    
    ##Prints the arc segment for the nozzle extension bell if div_r is specified by the user based on their bell geometry
    
    print(f'\t arc 14 15 ({aux_points[2,0]} {aux_points[2,1]} {aux_points[2,2]})')
    #edge5_text = captured_output.getvalue()



##Prints the arc segments at the throat section of the nozzle

edge3_2 = (f'\t arc 21 22 ({aux_points[3,0]} {aux_points[3,1]} {aux_points[3,2]})')
edge3_text.append(edge3_2)
edge3_text2 = "\n".join(edge3_text)
edge4_2 = (f'\t arc 23 24 ({aux_points[4,0]} {aux_points[4,1]} {aux_points[4,2]})')
edge4_text.append(edge4_2)
edge4_text2 = "\n".join(edge4_text)
#print(f'\t arc 25 26 ({aux_points[5,0]} {aux_points[5,1]} {aux_points[5,2]})')

print(");")



blockMeshDict = 'blockMeshDict.txt'


blockMeshDict_text = r'''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
  VERT_TEXT
);

blocks
(
    BLOX_TEXT
    hex (7 27 16 7 8 30 19 8 ) (20 1 800) simpleGrading (1 1 1)//
    hex (27 28 17 16 30 29 18 19 ) (30 1 800) simpleGrading (1 1 1)//
);

edges
(
    EDGE1_txt
    EDGE2_txt
    EDGE3_txt
    EDGE4_txt
    //EDGE5_txt
);

boundary
(
    inlet
    {
        type patch;
        faces
        (
            INLET_TEXT
            
        );
    }
    asym1
    {
        type wedge;
        faces
        (
            ASYM1_TEXT
            (7 16 19 8)//
            (16 17 18 19)//
        );
    }



    outlet
    {
        type patch;
        faces
        (
           //(7 27 16 7 ) 
            (8 30 19 8)//
            (30 29 18 19)//
            //(17 28 29 18) //
        );
    }

    outlet_r
    {
        type patch;
        faces
        (
           //(7 27 16 7 ) 
           // (8 30 19 8)//
            //(30 29 18 19)//
            (17 28 29 18) //
        );
    }



    asym2
    {
        type wedge;
        faces
        (
            ASYM2_TEXT
            (8 30 27 7)//
            (30 29 28 27) //
        );
    }
    nozzle
    {
        type wall;
        faces
        (
            NOZZLE_TEXT
            (16 27 28 17)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //'''



replacements = {"VERT_TEXT": vertices_text,
                "BLOX_TEXT": str(blocks_text2),
                "EDGE1_txt": str(edge12_text),
                "EDGE2_txt": str(edge2_text2),
                "EDGE3_txt": str(edge3_text2),
                "EDGE4_txt": str(edge4_text2),
                #"EDGE5_txt": str(edge5_text),
                "INLET_TEXT": str(inlet_text2),
                "ASYM1_TEXT": str(asym1_text2),
                "ASYM2_TEXT": str(asym2_text),
                "NOZZLE_TEXT": str(nozzle_text2)

}
for old, new in replacements.items(): ##this replaces the ambient pressure and chamber pressure specified by the user into the pressure file for OpenFOAM
    blockMeshDict_text = blockMeshDict_text.replace(old, new)

block_file = open(blockMeshDict, 'w')
block_file.write(blockMeshDict_text)
block_file.close()



folder_name = input("Name your OpenFOAM case: ")

##Creating the folder for the OpenFOAM case and placing all files in the correct directory for export
os.makedirs(folder_name, exist_ok=True)
os.makedirs(f"{folder_name}/0", exist_ok=True)
os.makedirs(f"{folder_name}/constant", exist_ok=True)
os.makedirs(f"{folder_name}/system", exist_ok=True)

##Move files before changing directory, removing .txt extension
for filename in ["U.txt", "T.txt", "alphat.txt", "epsilon.txt", "k.txt", "nut.txt", "p.txt"]:
    new_name = os.path.splitext(filename)[0]  # Removes .txt
    shutil.move(filename, os.path.join(folder_name, "0", new_name))

for filename in ["thermophysicalProperties.txt", "turbulenceProperties.txt"]:
    new_name = os.path.splitext(filename)[0]
    shutil.move(filename, os.path.join(folder_name, "constant", new_name))

for filename in ["controlDict.txt", "fvSchemes.txt", "fvSolution.txt","blockMeshDict.txt","decomposeParDict.txt"]:
    new_name = os.path.splitext(filename)[0]
    shutil.move(filename, os.path.join(folder_name, "system", new_name))

os.chdir(folder_name)
