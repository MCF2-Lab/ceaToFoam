#wrapper code for CEARUN using Python "CEA-Wrap" library using version 1.7.4 of CEA_Wrap
import shutil
import sys
from matplotlib import pyplot as plt
import numpy as np
from cantera import *
import cantera as ct
from CEAtoFOAM import Fuel, Oxidizer, RocketProblem, Mach_Exit, design_thrust,fuel_mdot,oxidizer_mdot ,Pe, Ve, mdot, A_star, Ae_At, Ae, Dt, De, Rt, Re, Ec, L_star, V_chamber, Ac, Dc, Rc, convergent_half_angle, divergent_half_angle, R_throat, Lconv, Ldiv, V_chamber_new, L_cylindrical,gamma,P1,p0,R,T1,OF_Ratio
import os
import math



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

comb_r=0.215 # radius of combustion chamber
comb_l=0.435 # straight length of combustion chamber
comb_r1=0.150 # curvature radius of combustion chamber tapering
comb_half_angle= 35 # half angle of the combustion chamber tapering
throat_r=0.0829 # radius at the nozzle
comb_full_l=0.723 # length of combustion chamber till throat
throat_curv_r=0.165 # curvature radius at the nozzle throat
div_half_angle_1=20 # half angle of the diverging section
straight_1=0.01 # first transitional straight section
div_r=5.8# curvature radius of the diverging nozzle bell
straight_2=0.014 # straight section at the nozzle exit
div_half_angle_2=10 # half angle of the diverging section at the termination
total_l=1.779 # total length of the engine
bell_r=0.360 # radius of the nozzle bell at the rim
domain_r=6*bell_r # radius of the area behind the nozzle for the plume
domain_ext=1.5*total_l+4 # extra length for the plume
wedge_half_angle=2.5 # half angle to create a wedge

points = np.zeros((31, 3)) # main blockmesh vertices
aux_points = np.zeros((6, 3)) # auxiliary points to specify arc segments

# points 0-8 are placed on the wedge axis

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
plt.xlim(left=-0.1)                # set up lower y-axis limit at zero
plt.xlim(right=10)
                   # set upper limit of the graph at 110% of maximum beam radius
plt.xlabel('Distance, [m]')
plt.ylabel('Radius, [m]')


#plt.savefig('Rothe_nozzle_profile.png', dpi=300)

# arc central point in combustion chamber
aux_points[0]=(comb_l+comb_r1*(math.sin(math.pi*comb_half_angle/360)),comb_r-comb_r1*(1-math.cos(math.pi*comb_half_angle/360)), 0)
#arc central point in the throat
aux_points[1]=(comb_full_l,throat_r, 0)

# arc central point in the diverging bell
midangle=(div_half_angle_1+div_half_angle_2)/2
x_center=points[14,0]+div_r*math.sin(math.pi*div_half_angle_1/180)
y_center=points[14,1]-div_r*math.cos(math.pi*div_half_angle_1/180)

aux_points[2]=(x_center-div_r*math.sin(math.pi*midangle/180), y_center+div_r*math.cos(math.pi*midangle/180 ),0)

#plt.scatter(aux_points[0:3,0], aux_points[0:3,1])
plt.show()
