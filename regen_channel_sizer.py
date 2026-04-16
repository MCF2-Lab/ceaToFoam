#based on Robert Watzlavick's regen engine sizer but makes sizing automatic based on the propellant choice, thrust, chamber pressure, and O/F ratio

from CEAtoFOAM import Ae, De, Dt, L_cylindrical, Lconv, OF_Ratio, R_throat, V_chamber_new, V_cone, gamma, chamber_temp, meanMolarMass, temp_fuel, temp_oxidizer, Ac_At, Ae_At, convergent_half_angle, divergent_half_angle, L_star, Dc, Ldiv, A_t, mdot,Mach_Exit


OFRatio = OF_Ratio #dimensionless
Ec = Ac_At #dimensionless (contraction ratio)
Ed = Ae_At #dimensionless (expansion ratio)
g = 32.2 #ft/s^2
Pc = float(input("Enter the chamber pressure in PSI: "))
Thrust = float(input("Enter the thrust in lbf: "))
Pe = float(input("Enter the exit pressure in PSI: "))
gammas = gamma #dimensionless
chamber_temperature = chamber_temp * 1.8 #Rankine
MolWeight = meanMolarMass #lb/mol
mdot = mdot*2.205 #kg/s to lbm/s

print(f"Calculated mean molar mass: {MolWeight:.2f} lb/mol")

L_star_print_out = """
            Propellant Combination              | L* (inches)
------------------------------------------------|----------------
Chlorine trifluoride/hydrazine-base fuel        | 20-35
Liquid fluorine/hydrazine                       | 24-28
Liquid fluorine/liquid hydrogen (GH2 injection) | 22-26
Liquid fluorine/liquid hydrogen (LH2 injection) | 25-30
Hydrogen peroxide/RP-1 (including catalyst bed) | 60-70
Nitric acid/hydrazine-base fuel                 | 30-35
Nitrogen tetroxide/hydrazine-base fuel          | 30-35
Liquid oxygen/ammonia                           | 30-40
Liquid oxygen/liquid hydrogen (GH2 injection)   | 22-28
Liquid oxygen/liquid hydrogen (LH2 injection)   | 30-40
Liquid oxygen/RP-1                              | 40-50
"""
print(L_star_print_out)


dpInj_fuel = 0.20*Pc #PSI

dpInj_oxidizer = (0.20*Pc)-2 #PSI

Cd = float(input("Enter the discharge coefficient (typically between 0.65 and 0.85): "))

tFuelAmbient = temp_fuel*1.8 #R
tOxAmbient = temp_oxidizer*1.8 #R

tFuelInjection = 100 + tFuelAmbient

tFuelAverage = (tFuelAmbient + tFuelInjection) / 2 #R

ox_density = 71.23 #lb/ft^3

R_universal = 1545.35 #ft-lbf/(lbmol-R) 
R_local = R_universal / meanMolarMass #ft-lbf/(lbmol-R)

fuelDensity_inj = tFuelInjection*(-0.02488)+63.63 #lb/ft^3
fuelDensity_ambient = tFuelAmbient*(-0.02488)+63.63 #lb/ft^3
fuelDensity_Average = tFuelAverage*(-0.02488)+63.63 #lb/ft^3

nPairs = float(input("Enter the number of injector pairs: "))

oxidizer_injection_angle = float(input("Enter the oxidizer injection angle in degrees (typically between 30 and 60 degrees): "))

#both the correction factors are based on the performance levels for a "well-designed" nozzle given by Huzel and Huang in their book "Modern Engineering for Design of Liquid-Propellant Rocket Engines".
C_star_correction_factor = 0.975
Cf_correction_factor = 0.980

angleConvNoz = convergent_half_angle #degrees
angleDivNoz = divergent_half_angle #degrees

L_star = L_star*39.37 #convert from meters to inches


chamber_diameter = Dc*39.37 #convert from meters to inches
throat_diameter = Dt*39.37 #convert from meters to inches
exit_diameter = De*39.37 #convert from meters to inches

chamber_length = L_cylindrical*39.37 #convert from meters to inches
convergent_cone_length = Lconv*39.37 #convert from meters to inches
divergent_cone_length = Ldiv*39.37 #convert from meters to inches

V_convergent = V_cone*61020 #convert from m^3 to in^3
V_chamber = V_chamber_new*61020 #convert from m^3 to in^3

Tc_ns = chamber_temperature

Tc_ns_design =  Tc_ns*((C_star_correction_factor)**2)

Tcoolant_init = tFuelAmbient

#allowable temperature of the hot wall will be 500F below the melting temperature of the material

Tmelting_material = float(input("Enter the melting temperature of the material in Rankine: ")) 
Tallowable_hot_wall = Tmelting_material - 500

stagnation_recovery_factor = 0.923 #also known as "R-factor" from Modern Engineering for Design of Liquid-Propellant Rocket Engines by Huzel and Huang

Btu_Conversion_Factor = 778.00 #ft-lbf/Btu,  J = energy conversion factor Pg 7 Huzel and Huang

Twc_ideal = Tallowable_hot_wall - 100

mu_gas = (4.66*(10**-9))*((meanMolarMass)**0.5)*(Tc_ns_design**0.6)

Prandtl_number_gas = (4*gamma)/((9*gamma)-5)

Cp_gas = (gamma*R_local)/(gamma-1) #specific heat of the gas mixture in the chamber

Cp_gas_BTU = Cp_gas / Btu_Conversion_Factor

Cp_fuel_avg = (tFuelAverage*0.00058009)+0.1598

Cp_fuel_inj = (tFuelInjection*0.00058009)+0.1598

Kin_Visc_fuel_avg = 0.0101187112-(0.0000897012025*tFuelAverage)+(0.000000342181242*(tFuelAverage**2))-(0.000000000725677457*(tFuelAverage**3))+(9.22297312E-13*(tFuelAverage**4))-(7.01648978E-16*(tFuelAverage**5))+(2.95605688E-19*(tFuelAverage**6))-(5.31734897E-23*(tFuelAverage**7))

Kin_Visc_fuel_inj = 0.0101187112-(0.0000897012025*tFuelInjection)+(0.000000342181242*(tFuelInjection**2))-(0.000000000725677457*(tFuelInjection**3))+(9.22297312E-13*(tFuelInjection**4))-(7.01648978E-16*(tFuelInjection**5))+(2.95605688E-19*(tFuelInjection**6))-(5.31734897E-23*(tFuelInjection**7))

Dyn_Visc_fuel_avg = Kin_Visc_fuel_avg * (fuelDensity_Average/12)

Dyn_Visc_fuel_inj = Kin_Visc_fuel_inj * (fuelDensity_inj/12)

k_fuel_avg = tFuelAverage*(-0.0000000013096)+0.000002233 #from Handbook of Aviation Fuel Properties

k_fuel_inj = tFuelInjection*(-0.0000000013096)+0.000002233 #from Handbook of Aviation Fuel Properties

k_AL = 220*(1.338*(10**-5)) #thermal conductivity of chamber wall material 

fin_factor = 0.600 #not sure where this came from but the citation given by Watzlavick was "Hill and Peterson pg 430". May be the Mechanics and Thermodynamics of Propulsion book? This is for the chamber wall material (in this )

R_throat = R_throat*39.37 #convert from meters to inches

hg_correction = 1.45 #from Watzlavick spreadsheet, its a correction factor based on past regen engine runs for the heat transfer coefficient


Chamber_temperature_1 = Tc_ns_design*stagnation_recovery_factor


number_tubes = float(input("Enter the number of cooling tubes in the regen channel: "))

M_e = Mach_Exit

M_chamber = 0.066

M_throat = 1.0

#breaking the denominator for the sigma correction factor into two pieces so that the terms can be evaluated separately in the main equation
sigma_chamber_sec1 = ((0.5*0.8*(1+(((gamma-1)/2)*(M_chamber**2))))+0.5)**0.68
sigma_chamber_sec2 = (1+(((gamma-1)/2)*(M_chamber**2)))**0.12
sigma_throat_sec1 = ((0.5*0.8*(1+(((gamma-1)/2)*(M_throat**2))))+0.5)**0.68
sigma_throat_sec2 = (1+(((gamma-1)/2)*(M_throat**2)))**0.12
sigma_exit_sec1 = ((0.5*0.8*(1+(((gamma-1)/2)*(M_e**2))))+0.5)**0.68
sigma_exit_sec2 = (1+(((gamma-1)/2)*(M_e**2)))**0.12



sigma_chamber = 1/(sigma_chamber_sec1*sigma_chamber_sec2) #where 0.8 is obtained from Fig 4-28 in Modern Engineering for Design of Liquid-Propellant Rocket Engines by Huzel and Huang
sigma_throat = 1/(sigma_throat_sec1*sigma_throat_sec2)
sigma_exit = 1/(sigma_exit_sec1*sigma_exit_sec2)

chamber_area = 3.14159*(chamber_diameter/2)**2
throat_area = 3.14159*(throat_diameter/2)**2
exit_area = 3.14159*(exit_diameter/2)**2

c_star = ((gamma*g*chamber_temperature*R_universal/MolWeight)**0.5)/(gamma*(((2/(gamma+1))**((gamma+1)/(gamma-1)))**0.5))*C_star_correction_factor #(Pc*throat_area)/mdot

hg_nominal = (0.026/(throat_diameter**0.2))*(((mu_gas**0.2)*Cp_gas_BTU)/((Prandtl_number_gas**0.6))*((Pc*g/c_star)**0.8)*((throat_diameter/R_throat)**0.1))


Ac_At = chamber_area/throat_area
Ae_At = exit_area/throat_area

Rd_Chamber = 1650.00
Rd_Throat = 1100.00
Rd_Exit = 1400.00

hg_chamber = 1/(1/(hg_nominal*(1/Ac_At)**0.9*sigma_chamber)+Rd_Chamber)*hg_correction #Btu/s-in2-R
hg_throat = 1/(1/(hg_nominal*sigma_throat)+Rd_Throat)*hg_correction #Btu/s-in2-R
hg_exit = 1/(1/(hg_nominal*(1/Ae_At)**0.9*sigma_exit)+Rd_Exit)*hg_correction #Btu/s-in2-R

T_chamber = Tc_ns_design*stagnation_recovery_factor

T_throat = T_chamber/(1+(gamma-1)/2*(M_throat**2))
T_exit = T_chamber/(1+(gamma-1)/2*(M_e**2))

print("---- HEAT TRANSFER CALCULATION RESULTS ----")

print("Mixture Ratio (O/F) = ", OFRatio)
print("Nominal Hot Wall Temperature = ", Twc_ideal, "R")
print("C* Value = ", c_star, "ft/s")
print("Chamber Pressure = ", Pc, "psi")
print("Exit Pressure = ", Pe, "psi")
print("Combustion Temperature (T0g) = ", chamber_temperature, "R")
print("Design Tc_ns (Stagnation Factor Included) = ", Tc_ns_design, "R")
print("Molecular Weight of the Combustion Products = ", MolWeight, "lb/mol")
print("L* = ", L_star, "in")
print("Injector Pressure Drop (Fuel) = ", dpInj_fuel, "psi")
print("Injector Pressure Drop (Oxidizer) = ", dpInj_oxidizer, "psi")
print("Discharge Coefficient = ", Cd)
print("Oxidizer Density = ", ox_density, "lb/ft^3")
print("Fuel Density = ", fuelDensity_Average, "lb/ft^3")
print("Number of Injector Pairs = ", nPairs)
print("Angle of Injector Orifice = ", oxidizer_injection_angle, "deg")
print("Specific Heat Ratio of Gases = ", gamma)
print("C* Correction Factor = ", C_star_correction_factor)
print("Thrust Coefficient Correction Factor = ", Cf_correction_factor)
print("Contraction Ratio = ", Ac_At)
print("Chamber Diameter = ", chamber_diameter, "in")
print("Throat Diameter = ", throat_diameter, "in")
print("Exit Diameter = ", exit_diameter, "in")
print("Chamber Length = ", chamber_length, "in")
print("Length of Convergent Section = ", convergent_cone_length, "in")
print("Length of Divergent Section = ", divergent_cone_length, "in")
print("Chamber Volume = ", V_chamber,"in^3")
print("Convergent Volume = ", V_convergent,"in^3")
print("Prandtl Number = ", Prandtl_number_gas)
print("Dynamic Viscosity of the Gas = ", mu_gas)
print("Specific Heat of the Gas at Constant Pressure = ", Cp_gas_BTU, "Btu/lbm-R")
print("Cp Fuel Avg = ", Cp_fuel_avg, "Btu/lbm-R")
print("Cp Fuel Injection = ", Cp_fuel_inj, "Btu/lbm-R")
print("Kinematic Viscosity of Fuel (Average) = ", Kin_Visc_fuel_avg, "ft^2/s")
print("Kinematic Viscosity of Fuel (Injection) = ", Kin_Visc_fuel_inj, "ft^2/s")
print("Thermal Conductivity of Fuel (Average)= ", k_fuel_avg, "Btu/s-in^2-R")
print("Thermal Conductivity of Fuel (Injection)= ", k_fuel_inj, "Btu/s-in^2-R")
print("Coolant Inlet Temperature = ", tFuelAmbient, "R")
print("Gas Constant of the Combustion Products = ", R_local, "ft-lbf/lbm-R")
print("Specific Heat of Gas Mixture = ", Cp_gas, "ft-lbf/lbm-R")
print("Thermal Conductivity of Material = ", k_AL, "Btu-in/s-in^2-R")
print("Fin Factor = ", fin_factor)
print("Radius of Curvature at Throat = ", R_throat, "in")
print("Number of Coolant Tubes = ", number_tubes)
print("Heat Transfer Coefficient Correction Factor = ", hg_correction)
print("Mach Number in Chamber = ", M_chamber)
print("Temperature in the chamber = ", T_chamber, "R")
print("Temperature at the throat = ", T_throat, "R")
print("Temperature at the exit = ", T_exit, "R")
print("Sigma Chamber = ", sigma_chamber)
print("Sigma Throat = ", sigma_throat)
print("Sigma Exit = ", sigma_exit)
print("Gas Side Convective Heat Transfer Coefficients = ", hg_chamber, hg_throat, hg_exit, "Btu/s-in2-R")


print("Allowable Hot Wall Temperature = ", Tallowable_hot_wall, "R")

print("Exit Mach Number = ", M_e)

print("Nominal gas side heat transfer coefficient = ", hg_nominal)

print("Gas side heat transfer coefficients", hg_chamber, hg_throat, hg_exit)

q_chamber = hg_chamber*(T_chamber-Tallowable_hot_wall)
q_throat = hg_throat*(T_throat-Tallowable_hot_wall)
q_exit = hg_exit*(T_exit-Tallowable_hot_wall)

print("q_chamber = ", q_chamber, "Btu/s-in2")
print("q_throat = ", q_throat, "Btu/s-in2")
print("q_exit = ", q_exit, "Btu/s-in2")

chamber_wall_thickness = float(input("Enter the wall thickness of the chamber in inches: "))
throat_wall_thickness = float(input("Enter the wall thickness of the throat in inches: "))
nozzle_wall_thickness = float(input("Enter the wall thickness of the nozzle in inches: "))

chamber_Twc = Tallowable_hot_wall-((q_chamber*chamber_wall_thickness)/k_AL)
throat_Twc = Tallowable_hot_wall-((q_throat*throat_wall_thickness)/k_AL)
nozzle_Twc = Tallowable_hot_wall-((q_exit*nozzle_wall_thickness)/k_AL)
print("Wall Temperature at the chamber = ", chamber_Twc, "R")
print("Wall Temperature at the throat = ", throat_Twc, "R")
print("Wall Temperature at the nozzle exit = ", nozzle_Twc, "R")

