#wrapper code for CEARUN using Python "CEA-Wrap" library using version 1.7.4 of CEA_Wrap
import shutil
import sys
from matplotlib import pyplot as plt
import numpy as np
from cantera import *
import cantera as ct
from CEA_Wrap import Fuel, Oxidizer, RocketProblem
import os
import re
from io import StringIO
from contextlib import redirect_stdout

plt.rcParams['axes.titlesize'] = 20  # Title font size
plt.rcParams['axes.labelsize'] = 20  # Axis label font size
plt.rcParams['xtick.labelsize'] = 16  # X-tick label font size
plt.rcParams['ytick.labelsize'] = 16  # Y-tick label font size
plt.rcParams['legend.fontsize'] = 14 # Legend font size
plt.rcParams['lines.linewidth'] = 1  # Line width for plots
plt.rcParams['lines.markersize'] = 8  # Marker size for points
plt.rcParams['figure.figsize'] = (10, 6)  # Default figure size
plt.rcParams['axes.grid'] = False  # Add grid to plots
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['legend.loc'] = 'upper left'  # Legend location


Pressure_Input = float(input("Enter Chamber Pressure (PSI): "))
Fuel_Input = input("Enter Fuel (RP-1,H2(L),CH4(L)): ")
Oxidizer_Input = input("Enter Oxidizer (O2(L),N2O4,O2): ")


fuel_temps = {
    "RP-1": 298.15,
    "H2(L)": 20.283,
    "CH4(L)": 111.66,
    "MMH": 292.0
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

OF_ratio = np.arange(1.0, 5.0 + 0.05, 0.05)  # start, stop (inclusive), step
Isp_values = []
temperatures = []
for OF in OF_ratio:
    mat1 = Fuel(str(Fuel_Input), temp_fuel, wt_percent=100, mols=None, chemical_composition = None, hf = None)
    mat2 = Oxidizer(str(Oxidizer_Input), temp_oxidizer, wt_percent=100, mols=None, chemical_composition=None, hf=None)
    problem = RocketProblem(pressure=Pressure_Input, massf=False, o_f=OF, pip = Pressure_Input/14.7)
    problem.run_cea(mat1, mat2)

    chamber_temp = None
    try:
        with open('my_output.out', 'r') as infile:
            for line in infile:
                if "T, K" in line:
                    nums = re.findall(r'[-+]?\d*\.\d+|\d+', line)
                    if nums:
                        # take last numeric token on the line (adjust if CEARUN format differs)
                        chamber_temp = float(nums[0])
                        print(f"Chamber Temperature for O/F={OF}: {chamber_temp} K")
                        break
    except FileNotFoundError:
        chamber_temp = float('nan')
        
    Isp = None
    try:
        with open('my_output.out', 'r') as infile:
            for line in infile:
                if "Isp" in line:
                    nums = re.findall(r'[-+]?\d*\.\d+|\d+', line)
                    if nums:
                        # take last numeric token on the line (adjust if CEARUN format differs)
                        Isp = float(nums[1])
                        print(f"Specific Impulse for O/F={OF}: {Isp} s")
                        break
    except FileNotFoundError:
        Isp = float('nan')

    # append once per OF iteration (use NaN if parsing failed)
    temperatures.append(chamber_temp)
    Isp_values.append(Isp)

# plot after the loop so x and y lengths match
plt.plot(OF_ratio, temperatures, marker='o', linestyle='-',color='black',markersize=3)
plt.plot(OF_ratio, Isp_values, marker='o', linestyle='--', color='black',markersize=3)
max_isp = np.max(Isp_values) if Isp_values else None
if max_isp is not None:
    print("Maximum Specific Impulse (s): ", max_isp)
    stoich_OF_Ratio_isp = OF_ratio[np.argmax(Isp_values)]
    print("Stoichiometric O/F Ratio for Max Isp: ", stoich_OF_Ratio_isp)
    plt.axvline(x=stoich_OF_Ratio_isp, color='black', linestyle='--', label='O/F Ratio for Max Isp')
max_temp = np.max(temperatures) #temperature of the stoichiometric mixture ratio (which is also identified through the graph)
print("Maximum Chamber Temperature (K): ", max_temp)
stoich_OF_Ratio = OF_ratio[np.argmax(temperatures)]
print("Stoichiometric O/F Ratio: ", stoich_OF_Ratio)
plt.axvline(x=stoich_OF_Ratio, color='black', linestyle=':', label='Stoichiometric O/F Ratio')


def on_pick(event):
    ind = event.ind[0] # Get the index of the picked point
    x_val = event.artist.get_xdata()[ind]
    y_val = event.artist.get_ydata()[ind]
    captured_output = StringIO()

# Use redirect_stdout to capture the output within a 'with' block
    with redirect_stdout(captured_output):
      print(f"User Picked O/F Ratio: ({x_val}, {y_val})")

# Get the captured output as a string
    selected_of_ratio = captured_output.getvalue()
    print(selected_of_ratio)  # Print the captured output to the console

x = OF_ratio
y = temperatures
line, = plt.plot(x, y, 'o', picker=5,color = 'black',markersize = 3) # Enable picking for the points

plt.gcf().canvas.mpl_connect('pick_event', on_pick)



plt.xlabel('O/F Ratio')
plt.ylabel('Chamber Temperature (K)')
plt.title(f'Chamber Temperature vs O/F Ratio for {Fuel_Input} and {Oxidizer_Input}')
plt.legend()
plt.show()
