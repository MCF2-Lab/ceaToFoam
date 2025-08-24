#wrapper code for CEARUN using Python "CEA-Wrap" library

from CEA_Wrap import Fuel, Oxidizer, RocketProblem

#using RP-1 and LOX as Propellants
#Chemical composition is not defined, thus, the program will use the default CEA values for the specified propellant

mat1 = Fuel("RP-1", temp=298.15, wt_percent=100, mols=None, chemical_composition = None, hf = None)
mat2 = Oxidizer("O2(L)",temp=90.170,wt_percent=100,mols=None,chemical_composition=None,hf=None)

#default unit for pressure is PSI
#massf is set to "True" so as to output the mass fractions of the reaction
#"pip" is the supersonic area ratio, which is the ratio of the exit area to the throat area, comprising the divergent section of the nozzle
problem = RocketProblem(pressure=300,massf=False,o_f=1.8,pip=20.4)
problem.run_cea(mat1,mat2)

#this ouputs the results of the CEARUN file for the LR-101 engine configuration
