#wrapper code for CEARUN using Python "CEA-Wrap" library using version 1.7.4 of CEA_Wrap
import sys
import numpy as np
from cantera import *
import cantera as ct
from CEA_Wrap import Fuel, Oxidizer, RocketProblem

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

import numpy as np

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

import numpy as np


def extract_m_inverse_n(input_file, output_file):
    """
    Extracts the three M(1/n) values from the output file, writes them to a text file,
    and returns them as a 3x1 list.
    """
    m_inverse_n_values = []
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.strip().startswith("M, (1/n)"): #locates the lines in the output file that denotes the molecular weight value in that section of the engine
                # Split the line and take the last three values
                parts = line.split()
                # Find the indices of the numeric values (skip 'M,' and '(1/n)')
                numeric_parts = [p for p in parts if p.replace('.', '', 1).isdigit()]
                # Convert to float and take the first three
                m_inverse_n_values = [float(numeric_parts[i]) for i in range(3)]
                break

    # Write to output file
    with open(output_file, 'w') as outfile:
        for value in m_inverse_n_values:
            outfile.write(f"{value}\n")

    return m_inverse_n_values

# Example usage:
m_inverse_n_list = extract_m_inverse_n('my_output.out', 'M_inverse_n.txt')
#print("Mixture Molecular Weight  [Chamber, Throat, Nozzle]:", m_inverse_n_list)
average_m_inverse_n = sum(m_inverse_n_list) / len(m_inverse_n_list)
print(round(float(average_m_inverse_n),2))  # computes and prints the average molecular weight in kg/kmol (g/mol = kg/kmol) according to the NASA CEA Analysis manual I


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

def format_q_from_file(input_file): #formats the mole fractions for each species involved in the reaction into the format that is inputted into the Sutherland coefficient section such that Cantera can read the information
    q_list = []  
    with open(input_file, 'r') as infile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) == 2:
                species, value = parts
                q_list.append(f"{species}:{value}")
    q_str = ' '.join(q_list)
    return f"q = {q_str}"
    #then formats the "q" string into a format such as this: q = H2:0.11111 O2:0.22222 H2O:0.33333 CO2:0.44444 for example


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
        gas.TPX = i, P1, q  # <-- FIXED
        cp1.append(gas.cp_mass)

    CpPoly1 = np.polyfit(Tv1, cp1, 4)
    
    c = tempRange2[0]
    d = tempRange2[1]

    Tv2 = np.arange(c,d,step)
    cp2 = []
        
    for i in range(c,d,step):
        gas.TPX = i, P1, q  # <-- FIXED
        cp2.append(gas.cp_mass)

    CpPoly2 = np.polyfit(Tv2, cp2, 4)

    return [Tv1, cp1, CpPoly1, Tv2, cp2, CpPoly2]

#######################################################################
# USER EDITABLE PART
#######################################################################
mech = 'gri30_highT.yaml'

P1 = Pressure_Input*6894.76; #convert PSI to Pa for use in the Sutherland coefficient calculation
T1 = 3538.04; #chamber temperature in Kelvin

Tlow = 200
Thigh = 6000
Tcommon = 1000
Tref = 298.15

q = q_formatted.replace("q = ", "").strip()  # Removes 'q = ' and any leading/trailing spaces
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

tempRange1 = [Tlow, Tcommon] #specifies the temperature range for the Openfoam file for the low and common temperature (which is usually a temperature somewhere in between the low and high value )
tempRange2 = [Tcommon, Thigh] #specifies the common temperature and the high temperature for the openfoam file. For the high temperature, it is likely the best practice to set it +3000K higher than the actual chamber temperature outputted by CEA.
step = 1
out = cpPolynomials(P1,q,mech,tempRange1,tempRange2, step)

R = GasConstant/(round(meanMolarMass,2)) #finds the specific gas constant for the mixture

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

#openfoam formatting from the output of the Cantera calculations
print("    specie")
print("    {")
print("        nMoles          1;")
print("        molWeight       %s;" % meanMolarMass)
print("    }")
print("    thermodynamics")
print("    {")
print("        Tlow            %f;" % Tlow)
print("        Thigh           %f;" % Thigh)
print("        Tcommon         %f;" % Tcommon)
print("        highCpCoeffs    ( %g %g %g %g %g %g %g );" % (Hcof[0],  Hcof[1], Hcof[2], Hcof[3], Hcof[4], high_enthalpy_offset,  high_entropy_offset))
print("        lowCpCoeffs     ( %g %g %g %g %g  %g %g );" % (Lcof[0], Lcof[1], Lcof[2], Lcof[3], Lcof[4],  low_enthalpy_offset, low_entropy_offset))
print("    }")
print("    transport")
print("    {")
print("        As              1.67212e-06;")
print("        Ts              170.672;")
print("    }")
print("}")

