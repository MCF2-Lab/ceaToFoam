#wrapper code for CEARUN using Python "CEA-Wrap" library using version 1.7.4 of CEA_Wrap
import shutil
import sys
import numpy as np
from cantera import *
import cantera as ct
from CEA_Wrap import Fuel, Oxidizer, RocketProblem
import os

#using RP-1 and LOX as Propellants
#Chemical composition is not defined, thus, the program will use the default CEA values for the specified propellant

Pressure_Input = float(input("Enter Chamber Pressure (PSI): "))
Fuel_Input = input("Enter Fuel (RP-1,H2(L),CH4(L)): ")
Oxidizer_Input = input("Enter Oxidizer (O2(L),N2O4,O2): ")
OF_Ratio = float(input("Enter O/F Ratio: "))

fuel_temps = {
    "RP-1": 298.15,
    "H2(L)": 20.283,
    "CH4(L)": 111.66
}

oxidizer_temps = {
    "O2(L)": 90.170,
    "N2O4": 294.15,
    "O2": 298.15
}

temp_fuel = fuel_temps.get(Fuel_Input)
temp_oxidizer = oxidizer_temps.get(Oxidizer_Input)

if temp_fuel is None or temp_oxidizer is None:
    print("Invalid fuel or oxidizer selection.")
    sys.exit(1)

mat1 = Fuel(str(Fuel_Input), temp_fuel, wt_percent=100, mols=None, chemical_composition = None, hf = None)
mat2 = Oxidizer(str(Oxidizer_Input),temp_oxidizer,wt_percent=100,mols=None,chemical_composition=None,hf=None)

#default unit for pressure is PSI
#massf is set to "True" so as to output the mass fractions of the reaction
#"pip" is the supersonic area ratio, which is the ratio of the exit area to the throat area, comprising the divergent section of the nozzle
problem = RocketProblem(pressure=Pressure_Input, massf=False, o_f=OF_Ratio, pip=18.9)
problem.run_cea(mat1, mat2)

#this ouputs the results of the CEARUN file for the LR-101 engine configuration
#my_output.inp is the input file containing the information entered above by the user
#my_output.out is the output file containing the results of the CEA run
#my_output.plt outputs the pressure, temperature, and other properties at the chamber, throat, and nozzle

def extract_rows_in_range(input_file, output_file, start_row, end_row):
    #if you want to run a new propellant combination, you must change the inputs of within the CEA_Code.py file to
    #the specifications of your engine, this file will handle the output from my_output.out and parse it into a format
    #which is applicable to the OpenFOAM simulation
    """
    Extracts rows within a specified range from the input file and writes them to the output file.

    :param input_file: Path to the input text file.
    :param output_file: Path to the output text file.
    :param start_row: Starting row number (1-based index).
    :param end_row: Ending row number (1-based index).
    """
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
    
    # Extract the specified range of rows
    extracted_lines = lines[start_row - 1:end_row]
    
    with open(output_file, 'w') as outfile:
        outfile.writelines(extracted_lines)

#this function extracts the chamber temperature from the output file
def extract_chamber_temperature(input_file):
    """
    Extracts the chamber temperature (first value) from the line containing 'T, K'.
    Returns it as a float.
    """
    with open(input_file, 'r') as infile:
        for line in infile:
            if "T, K" in line:
                parts = line.split()
                for part in parts:
                    try:
                        chamber_temp = float(part)
                        return chamber_temp
                    except ValueError:
                        continue
    #returns nothing if the chamber temperature is not found
    return None 

chamber_temp = extract_chamber_temperature('my_output.out')
print("Chamber temperature:", chamber_temp)


def extract_gammas(input_file):
    """
    Extracts the chamber temperature (first value) from the line containing 'GAMMAs'.
    Returns it as a float.
    """
    with open(input_file, 'r') as infile:
        for line in infile:
            if "GAMMAs" in line:
                parts = line.split()
                for part in parts:
                    try:
                        gamma = float(part)
                        return gamma
                    except ValueError:
                        continue
    #returns nothing if the chamber temperature is not found
    return None 

chamber_temp = extract_chamber_temperature('my_output.out')
print("Chamber temperature:", chamber_temp)

gammas = extract_gammas('my_output.out')
print("Gammas:", gammas)

def extract_m_inverse_n(input_file, output_file):
    """
    Extracts the three M(1/n) values from the output file, writes them to a text file,
    and returns them as a 3x1 list.
    """
    m_inverse_n_values = []
    with open(input_file, 'r') as infile:
        for line in infile:
            #locates the lines in the output file that denotes the molecular weight value in that section of the engine
            if line.strip().startswith("M, (1/n)"):
                # Split the line and take the last three values
                parts = line.split()
                # Finds the indices of the numeric values (skip 'M,' and '(1/n)')
                numeric_parts = [p for p in parts if p.replace('.', '', 1).isdigit()]
                # Convert to float and take the first three
                m_inverse_n_values = [float(numeric_parts[i]) for i in range(3)]
                break

    # Writes to an output file
    with open(output_file, 'w') as outfile:
        for value in m_inverse_n_values:
            outfile.write(f"{value}\n")

    return m_inverse_n_values


m_inverse_n_list = extract_m_inverse_n('my_output.out', 'M_inverse_n.txt')

#print("Mixture Molecular Weight  [Chamber, Throat, Nozzle]:", m_inverse_n_list)
average_m_inverse_n = sum(m_inverse_n_list) / len(m_inverse_n_list)

# computes and prints the average molecular weight in kg/kmol (g/mol = kg/kmol) according to NASA CEA Analysis manual I
print(f"Molecular Weight: {round(float(average_m_inverse_n),2)}")  


# -----------------------------------------------------------------------#
# The following code generates the thermodynamic coefficients for the OpenFOAM simulation
# using Cantera to compute the NASA polynomials for the given mixture at the chamber conditions
# -----------------------------------------------------------------------#
GasConstant = 8314.4621 # J/(mol*K)
meanMolarMass = round(average_m_inverse_n,2)

def extract_species_and_chamber_mole_fractions(input_file, output_file): 
  
    collecting = False
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not collecting and "MOLE FRACTIONS" in line.upper():
                collecting = True
                continue
            if collecting:
                if "THERMODYNAMIC PROPERTIES FITTED TO 20000.K" in line.upper():
                    break
                parts = line.strip().split()
                if len(parts) >= 2:
                    species = parts[0].replace('*', '')
                    try:
                        mole_fraction = float(parts[1])
                        outfile.write(f"{species} {mole_fraction:.5f}\n")
                    except ValueError:
                        continue
#formats the mole fractions for each species involved in the reaction into the format that is inputted into the Sutherland coefficient section such that Cantera can read the information
def format_q_from_file(input_file):
    q_list = []  
    with open(input_file, 'r') as infile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) == 2:
                species, value = parts
                q_list.append(f"{species}:{value}")
    q_str = ' '.join(q_list)
     #then formats the "q" string into a format such as this: q = H2:0.11111 O2:0.22222 H2O:0.33333 CO2:0.44444 for example
    return f"q = {q_str}"

extract_species_and_chamber_mole_fractions('my_output.out', 'mole_fractions_section.txt')
extract_species_and_chamber_mole_fractions('my_output.out', 'chamber_mole_fractions_only.txt')

q_formatted = format_q_from_file('chamber_mole_fractions_only.txt')

meanMolarMass = round(average_m_inverse_n,2) #g/mol = kg/kmol

def cpPolynomials(P1,q_new,mech,tempRange1,tempRange2,step):

    gas = ct.Solution(mech)
    gas.TPX = T1, P1, q_new  # <-- FIXED

    a = tempRange1[0]
    b = tempRange1[1]

    Tv1 = np.arange(a,b,step)
    cp1 = []
        
    for i in range(a,b,step):
        #these are fixed values
        gas.TPX = i, P1, q
        cp1.append(gas.cp_mass)

    CpPoly1 = np.polyfit(Tv1, cp1, 4)
    
    c = tempRange2[0]
    d = tempRange2[1]

    Tv2 = np.arange(c,d,step)
    cp2 = []
        
    for i in range(c,d,step):
        gas.TPX = i, P1, q
        cp2.append(gas.cp_mass)

    CpPoly2 = np.polyfit(Tv2, cp2, 4)

    return [Tv1, cp1, CpPoly1, Tv2, cp2, CpPoly2]

#######################################################################
# USER EDITABLE PART
#######################################################################
mech = 'gri30_highT.yaml'

 #convert PSI to Pa for use in the Sutherland coefficient calculation
P1 = Pressure_Input*6894.76;
 #chamber temperature in Kelvin
T1 = chamber_temp
gamma = gammas

#these are the Tlow and Thigh values for the OpenFOAM thermophysicalProperties file also coinciding with the Sutherland gri30_highT.yaml file temperature value ranges
Tlow = 200
Thigh = 6000 
Tcommon = 1000
Tref = 298.15

# Removes 'q = ' and any leading/trailing spaces
q = q_formatted.replace("q = ", "").strip()  
print(q)

gas = ct.Solution('gri30_highT.yaml')
allowed_species = set(gas.species_names)

def format_q_from_file_filtered(input_file, allowed_species):
    q_list = []
    with open(input_file, 'r') as infile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) == 2:
                species, value = parts
                if species in allowed_species:
                    q_list.append(f"{species}:{value}")
    q_str = ' '.join(q_list)
    return q_str

# Use the filtered function
q = format_q_from_file_filtered('chamber_mole_fractions_only.txt', allowed_species)

# Now q only contains valid species for Cantera
gas.TPX = T1, P1, q

#######################################################################
# END OF USER EDITABLE PART
#######################################################################


gas = ct.Solution(mech) 
gas.TPX = T1, P1, q  # <-- FIXED

#specifies the temperature range for the Openfoam file for the low and common temperature (which is usually a temperature somewhere in between the low and high value )
tempRange1 = [Tlow, Tcommon] 

 #specifies the common temperature and the high temperature for the openfoam file. For the high temperature, it is likely the best practice to set it +3000K higher than the actual chamber temperature outputted by CEA.
tempRange2 = [Tcommon, Thigh] 
step = 1
out = cpPolynomials(P1,q,mech,tempRange1,tempRange2, step)

#finds the specific gas constant for the mixture
R = GasConstant/(round(meanMolarMass,2)) 

Lcof_rev = out[2]/R
Hcof_rev = out[5]/R

Lcof = list(reversed(Lcof_rev))
Hcof = list(reversed(Hcof_rev))


hLoff = Lcof[0] + Lcof[1]*((Tref**1)/2) + Lcof[2]*((Tref**2)/3) + Lcof[3]*((Tref**3)/4) + Lcof[4]*((Tref**4)/5)
sLoff = Lcof[0]*np.log(Tref) + Lcof[1]*((Tref**1)/1) + Lcof[2]*((Tref**2)/2) + Lcof[3]*((Tref**3)/3) + Lcof[4]*((Tref**4)/4)

hHoff = Hcof[0] + Hcof[1]*((Tref**1)/2) + Hcof[2]*((Tref**2)/3) + Hcof[3]*((Tref**3)/4) + Hcof[4]*((Tref**4)/5)
sHoff = Hcof[0]*np.log(Tref)  + Hcof[1]*((Tref**1)/1) + Hcof[2]*((Tref**2)/2) + Hcof[3]*((Tref**3)/3) + Hcof[4]*((Tref**4)/4)

low_enthalpy_offset = gas.enthalpy_mass/R - hLoff*Tref
low_entropy_offset = gas.entropy_mass/R - sLoff

high_enthalpy_offset = gas.enthalpy_mass/R - hHoff*Tref
high_entropy_offset = gas.entropy_mass/R - sHoff

#the following section is the generation of all necessary files required for an OpenFOAM simulation of your engine according to the prescribed thermophysical properties previously inputted

a_sound = np.sqrt(gas.cp_mass*R* T1)

#Asssuming a Mach number of 0.1 in the chamber for the velocity input and finding the velocity at the inlet by multiplying the speed of sound by this Mach.
VELOCITY_INPUT = 0.1*a_sound
U = "U.txt"
U_text = '''/*--------------------------------*- C++ -*----------------------------------*\
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

internalField   uniform (XVEL 0 0);

boundaryField
{
    inlet
    {
        type            zeroGradient;
    }

    freestream
    {
        type            zeroGradient;
    }

    walls
    {
        type            noSlip;
    }

    bottomEmptyFaces
    {
        type            empty;
    }

    topEmptyFaces
    {
        type            empty;
    }
}


// ************************************************************************* //'''

U_text = U_text.replace("XVEL", VELOCITY_INPUT) #replaces the XVEL in the U file with the user inputted value for the velocity in the x-direction

U_file = open(U, 'w')
U_file.write(U_text)
U_file.close()


TEMPERATURE_INPUT = str(T1)
T = "T.txt"
GAMMA_INPUT = str(gamma)
T_text = '''/*--------------------------------*- C++ -*----------------------------------*\
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

internalField   uniform 300;
boundaryField
{
    outlet
    {
        type            waveTransmissive;
        gamma           GAMMA_VAL;
        fieldInf        300;
        lInf            10;
        value           uniform 300;
    }
    inlet
    {
        type         totalTemperature;
	    gamma        GAMMA_VAL;
        T0           uniform 3538.04;
	    value	     uniform 3538.04;
    }
    wall
    {
        type            zeroGradient;
    }
    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //'''

T_text = T_text.replace("T_VAL", TEMPERATURE_INPUT) #replaces the TVAL in the T file with the user inputted value for the temperature
T_text = T_text.replace("GAMMA_VAL", GAMMA_INPUT) #replaces the GAMMA_VAL in the T file with the user inputted value for the gamma

T_file = open(T, 'w')
T_file.write(T_text)
T_file.close()





alphat = "alphat.txt"
alphat_text = '''
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

    wall
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }

    frontAndBackPlanes
    {
        type            empty;
    }

}


// ************************************************************************* //'''


alphat_file = open(alphat, 'w')
alphat_file.write(alphat_text)
alphat_file.close()





epsilon = "epsilon.txt"
epsilon_text = '''
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

    wall
    {
        type            epsilonWallFunction;
        value           uniform 266000;
    }

    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //'''

epsilon_file = open(epsilon, 'w')
epsilon_file.write(epsilon_text)
epsilon_file.close()




k = "k.txt"

k_text = '''
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

    wall
    {
        type            kqRWallFunction;
        value           uniform 1000;
    }

    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //'''

k_file = open(k, 'w')
k_file.write(k_text)
k_file.close()



nut = "nut.txt"

nut_text = '''
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

    wall
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //'''

nut_file = open(nut, 'w')
nut_file.write(nut_text)
nut_file.close()


pressure = "p.txt"
ambient_pressure = input("Enter Ambient Pressure (PSI): ")
ambient_pressure_Pa = str(float(ambient_pressure) * 6894.76)  # Convert PSI to Pa

pressure_text = '''
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
    wall
    {
        type            zeroGradient;
    }
    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //'''

replacements = {"AMB_PRESSURE": ambient_pressure_Pa,
                "P_VAL": str(P1)
                ,"GAMMA_VAL": GAMMA_INPUT
}
for old, new in replacements.items(): #this replaces the ambient pressure and chamber pressure specified by the user into the pressure file for OpenFOAM
    pressure_text = pressure_text.replace(old, new)

pressure_file = open(pressure, 'w')
pressure_file.write(pressure_text)
pressure_file.close()




thermophysicalProperties = "thermophysicalProperties.txt"

tpp_text = '''/*--------------------------------*- C++ -*----------------------------------*\
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





for old, new in replacements_tpp.items(): #this replaces the ambient pressure and chamber pressure specified by the user into the pressure file for OpenFOAM
    tpp_text = tpp_text.replace(old, new)

thermo_file = open(thermophysicalProperties, 'w')
thermo_file.write(tpp_text)
thermo_file.close()


turbulenceProperties = "turbulenceProperties.txt"

turbProp_text = '''
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

control_dict_text = '''
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

fv_schemes_text = '''
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

fv_solutions_text = '''
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


folder_name = input("Name your OpenFOAM case: ")

# Creating the folder for the OpenFOAM case and placing all files in the correct directory for export
os.makedirs(folder_name, exist_ok=True)
os.makedirs(f"{folder_name}/0", exist_ok=True)
os.makedirs(f"{folder_name}/constant", exist_ok=True)
os.makedirs(f"{folder_name}/system", exist_ok=True)

# Move files before changing directory, removing .txt extension
for filename in ["U.txt", "T.txt", "alphat.txt", "epsilon.txt", "k.txt", "nut.txt", "p.txt"]:
    new_name = os.path.splitext(filename)[0]  # Removes .txt
    shutil.move(filename, os.path.join(folder_name, "0", new_name))

for filename in ["thermophysicalProperties.txt", "turbulenceProperties.txt"]:
    new_name = os.path.splitext(filename)[0]
    shutil.move(filename, os.path.join(folder_name, "constant", new_name))

for filename in ["controlDict.txt", "fvSchemes.txt", "fvSolution.txt"]:
    new_name = os.path.splitext(filename)[0]
    shutil.move(filename, os.path.join(folder_name, "system", new_name))

os.chdir(folder_name)
