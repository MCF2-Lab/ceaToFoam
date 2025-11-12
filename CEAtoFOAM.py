#wrapper code for CEARUN using Python "CEA-Wrap" library using version 1.7.4 of CEA_Wrap
import shutil
import sys
from matplotlib import pyplot as plt
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
Design_altitude = float(input("Enter Design Altitude (ft): "))
altitude = Design_altitude*0.3048 #converts feet to meters




#AMBIENT PRESSURE CALCULATION SECTION

#----------------------------------------------------------------------------------------#
if (11000*0.3048<altitude) and (altitude<25000*0.3048):
    Tamb = -56.46 #CELSIUS
    p0 = 1000*(22.65*np.exp(1.73-0.000157*altitude)) #Pascals, converting from kPa to Pa with *1000 https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
elif altitude>=25000*0.3048:
    Tamb = -131.21 + 0.00299*altitude #CELSIUS
    p0 = 1000*(2.488*((Tamb+273.1)/216.6)**(-11.388)) # the altitude pressure in Pascals
else:
    #For altitudes 0 to 11,000 m
    Tamb = 15.04 - 0.00649*altitude
    p0 = 1000*(101.29*((Tamb+273.1)/288.08)**5.256)
#---------------------------------------------------------------------------------------#


p0_psi = p0*0.000145038 #converts Pascals to PSI for use in the CEA code
ambient_pressure = p0_psi

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
problem = RocketProblem(pressure=Pressure_Input, massf=False, o_f=OF_Ratio, pip = Pressure_Input/ambient_pressure)
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

#these are the Tlow and Thigh values for the Openfoam thermophysicalProperties file also coinciding with the Sutherland gri30_highT.yaml file temperature value ranges
Tlow = 200
Thigh = 6000 
Tcommon = 1000
Tref = 298.15

# Removes 'q = ' and any leading/trailing spaces
q = q_formatted.replace("q = ", "").strip()  

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

print(q)
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

#this plots the Cp(T) vs Static Temperature (K) curve for the given mixture
npoints = 50
temRL = np.linspace(200, 1200, npoints)
temRH = np.linspace(800, 5200, npoints)

cpL = np.zeros(npoints)
cpH = np.zeros(npoints)

cpL = np.polyval(out[2], temRL)
cpH = np.polyval(out[5], temRH)

#For OpenFoam
cpLbyR = np.zeros(npoints)
cpHbyR = np.zeros(npoints)

cpLbyR = np.polyval(out[2], temRL)
cpHbyR = np.polyval(out[5], temRH)

plt.plot(temRL,cpLbyR, lw=2,color = 'black')
plt.plot(temRH,cpHbyR, lw=2,color = 'black',linestyle='--')

plt.xlabel('Temperature (K)')
plt.ylabel('$C_p$ [J/kg/K]')
plt.grid(color='black', alpha=0.5, linewidth=0.5)

#plt.savefig("cp_vs_T.png")
plt.show()

cpL = np.polyval(out[2], temRL)
cpH = np.polyval(out[5], temRH)

# specific gas constant used earlier
# R = GasConstant/(round(meanMolarMass,2))

# Cp/R (dimensionless-ish) for the same temperature ranges
cpL_divR = cpL / R
cpH_divR = cpH / R

# Plot Cp and Cp/R side-by-side
plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.plot(temRL, cpL, lw=2, label='Cp (low range)',  color = 'black')
plt.plot(temRH, cpH, lw=2, label='Cp (high range)', color = 'black',linestyle='--')
plt.xlabel('Temperature (K)')
plt.ylabel('Cp [J/kg/K]')
plt.title('Cp vs T')
plt.grid(alpha=0.4)
plt.legend()

plt.subplot(1,2,2)
plt.plot(temRL, cpL_divR, lw=2, label='Cp/R (low range)', color = 'black')
plt.plot(temRH, cpH_divR, lw=2, label='Cp/R (high range)',color = 'black', linestyle='--')
plt.xlabel('Temperature (K)')
plt.ylabel('Cp / R')
plt.title('Cp/R vs T')
plt.grid(alpha=0.4)
plt.legend()

plt.tight_layout()
plt.savefig("cp_and_cpOverR_vs_T.png", dpi=200, bbox_inches='tight')
plt.show()



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

a_sound = np.sqrt(gammas*R*T1)
print("speed of sound of the gas mixture in the combustion chamber:", a_sound)



#CHAMBER AND NOZZLE DIMENSION CALCULATIONS (Different from CEA, uses Isentropic flow equations and thus assumes a constant gamma starting from the chamber combustion point
#listed in the CEA output file, all the way to the nozzle exit plane point. This is a common simplification made by amateur rocket engine designers to get a good estimation of rocket
#engine chamber and nozzle dimensions without having to do complex thermodyncamic calculations at each point along the nozzle length, which
#takes shifting equilibrium chemistry into account as well as variable gamma values along the nozzle length
 

Mach_Exit = np.sqrt((2.0 / (gamma - 1.0)) * ((P1 / p0) ** ((gamma - 1.0) / gamma) - 1.0))

design_thrust = float(input("Enter desired thrust (lbs): "))
design_thrust = design_thrust * 4.44822  #converts thrust from lbs to Newtons

Pe = P1 * (1.0 + ((gamma - 1.0) / 2.0) * Mach_Exit ** 2) ** (-gamma / (gamma - 1.0))

Ve = np.sqrt((2.0 * gamma * R * T1) / (gamma - 1.0) * (1.0 - (Pe / P1) ** ((gamma - 1.0) / gamma)))

mdot = design_thrust / Ve  #kg/s

fuel_mdot = mdot/(1+OF_Ratio)
oxidizer_mdot = mdot-fuel_mdot

print(f"Mass Flow Rate (kg/s): {float(mdot)}")
# Use P1 (Pa) for numeric calculations and p0 (Pa) for ambient pressure
A_star = A_t = (mdot / P1) * np.sqrt((R * T1)/gamma) * ((((gamma-1)/2 + 1) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))))
 
 
Ae_At = (1.0 / Mach_Exit) * ((2.0 / (gamma + 1.0)) * (1.0 + ((gamma - 1.0) / 2.0) * Mach_Exit ** 2)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
 


Ae = Ae_At * A_star
 
Pe = P1 * (1.0 + ((gamma - 1.0) / 2.0) * Mach_Exit ** 2) ** (-gamma / (gamma - 1.0))
 
Ve = np.sqrt((2.0 * gamma * R * T1) / (gamma - 1.0) * (1.0 - (Pe / P1) ** ((gamma - 1.0) / gamma)))
 
Dt = 2*np.sqrt(A_t/np.pi)
De = 2*np.sqrt(Ae/np.pi)

Rt = Dt/2
Re = De/2

Dt_cm = Dt*100

Ec = Ac_At = (8*(Dt_cm)**-0.6)+1.25

L_star = 1.2
V_chamber = L_star * A_t

Ac = Ac_At * A_t

Dc = 2*np.sqrt(Ac/np.pi)
Rc = Dc/2

convergent_half_angle = float(input("Enter convergent half angle (degrees): "))
convergent_half_angle = np.radians(convergent_half_angle)
divergent_half_angle = float(input("Enter divergent half angle (degrees): "))
divergent_half_angle = np.radians(divergent_half_angle)

R_throat = 1.5*Rt

Lconv = (Rt*(np.sqrt(Ac_At)-1)+(R_throat)*(1/np.cos(convergent_half_angle)-1))/np.tan(convergent_half_angle)

Ldiv = (Rt*(np.sqrt(Ae_At)-1)+(R_throat)*(1/np.cos(divergent_half_angle)-1))/np.tan(divergent_half_angle)

V_cone = (1/3)*np.pi*(Rc**2+(Rc*Rt)+Rt**2)*Lconv

V_chamber_new = V_chamber - V_cone

L_cylindrical = V_chamber_new/Ac

total_engine_length = L_cylindrical + Lconv + Ldiv

#Output of the percentage error between the farfield ambient pressure and the exit pressure of the nozzle
Percent_Error_P_Exit = abs((Pe - p0)/p0)*100
print(f"Percent Error in Exit Pressure (%): {float(Percent_Error_P_Exit)}") 

print(f"Throat Area (m^2): {float(A_star)}")
print(f"Exit Area (m^2): {float(Ae)}")
print(f"Chamber Area (m^2): {float(Ac)}")
print(f"Throat Diameter (m): {float(Dt)}")
print(f"Throat Radius (m): {float(Rt)}")
print(f"Chamber Diameter (m): {float(Dc)}")
print(f"Chamber Radius (m): {float(Rc)}")
print(f"Exit Diameter (m): {float(De)}")
print(f"Exit Radius (m): {float(Re)}")
print("Contraction Ratio (Ec): {:.2f}".format(Ec))
print("Expansion Ratio (Ae/At): {:.2f}".format(Ae_At))
print(f"Exit Pressure (Pa): {float(Pe)}")
print(f"Exit Velocity (m/s): {float(Ve)}")
print(f"Exit Mach Number: {float(Mach_Exit)}")
print(f"Area Ratio (Ae/At): {float(Ae_At)}")
print(f"Mass Flow Rate (kg/s): {float(mdot)}")
print(f"Fuel Mass Flow Rate (kg/s): {float(fuel_mdot)}")
print(f"Oxidizer Mass Flow Rate (kg/s): {float(oxidizer_mdot)}")
print(f"Chamber Diameter (m): {float(Dc)}")
print(f"Convergent Length (m): {float(Lconv)}")
print(f"Cylindrical Chamber Length (m): {float(L_cylindrical)}")
print(f"Divergent Length (m): {float(Ldiv)}")
print("Chamber Volume (m^3): {:.4f}".format(V_chamber_new))
print(f"Radius of Curvature at Throat (m): {float(R_throat)}")
print(f"Chamber Temperature (K): {float(T1)}")
print(f"Total Engine Length (m): {float(total_engine_length)}")


def plot_nozzle_contour_piecewise(savefile='nozzle_contour.png'):
    """
    Piecewise-linear nozzle contour:
      - straight cylindrical chamber from x=0 to x=L_cylindrical (radius = Rc)
      - straight convergent from x=L_cylindrical to throat at x=L_cylindrical+Lconv (radius = Rt)
      - optional tiny throat smoothing (disabled here)
      - straight divergent cone from throat to exit (x=L_cylindrical+Lconv+Ldiv, radius=Re)
    This ensures the chamber section up to the first blue point is flat, and the diverging
    section is a straight line up to the exit point.
    """
    # key axial positions
    x0 = 0.0
    x_ch_end = L_cylindrical
    x_throat = L_cylindrical + Lconv
    x_exit = L_cylindrical + Lconv + Ldiv
    y_farfield = Re*1.5  # for plotting only
    x_farfield = 5*L_star  # arbitrary farfield location

    # piecewise sampling
    n_ch = 4
    n_conv = 80
    n_div = 80

    x_ch = np.linspace(x0, x_ch_end, n_ch)
    y_ch = np.full_like(x_ch, Rc)

    x_conv = np.linspace(x_ch_end, x_throat, n_conv)
    y_conv = np.linspace(Rc, Rt, n_conv)

    x_div = np.linspace(x_throat, x_exit, n_div)
    y_div = np.linspace(Rt, Re, n_div)

    # combine
    x = np.concatenate([x_ch, x_conv[1:], x_div[1:]])
    y = np.concatenate([y_ch, y_conv[1:], y_div[1:]])

    plt.figure(figsize=(9,4.5))
    plt.plot(x, y, '-r', lw=2, label='nozzle contour (piecewise linear)')
    plt.plot(x, -y, '-r', lw=2)  # mirrored lower half if desired
    # mark the key blue points: chamber end, throat, exit
    plt.scatter([x_ch_end, x_throat, x_exit], [Rc, Rt, Re], c='b', zorder=10)
    plt.annotate('chamber end', (x_ch_end, Rc), xytext=(6,6), textcoords='offset points')
    plt.annotate('throat', (x_throat, Rt), xytext=(6,-12), textcoords='offset points')
    plt.annotate('exit', (x_exit, Re), xytext=(6,6), textcoords='offset points')

    plt.xlabel('Axial distance [m]')
    plt.ylabel('Radius [m]')
    plt.title('Nozzle contour (chamber - convergent - throat - diverging cone)')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(savefile, dpi=200, bbox_inches='tight')
    plt.show()
    print(f'Nozzle contour plotted and saved to: {savefile}')

# replace call to previous plot function with this one
plot_nozzle_contour_piecewise()





 
