import os
import shutil
import numpy as np
import math
import matplotlib                       # import matplotlib and its its components for plotting
import matplotlib.pyplot as plt
from CEAtoFOAM import Rt,Rc,L_cylindrical,Ldiv,Lconv,R_throat,divergent_half_angle,convergent_half_angle,Re,folder_name
import io



comb_r=Rc # radius of combustion chamber
comb_l=L_cylindrical # straight length of combustion chamber
comb_r1=0.01# curvature radius of combustion chamber tapering
comb_half_angle= convergent_half_angle # half angle of the combustion chamber tapering
throat_r=Rt # radius at the nozzle
comb_full_l=L_cylindrical+Lconv # length of combustion chamber till throat
throat_curv_r=R_throat # curvature radius at the nozzle throat
div_half_angle_1=divergent_half_angle # half angle of the diverging section
straight_1=0.01 # first transitional straight section
div_r = None  # no bell curvature for a purely conical nozzle
straight_2=0.014 # straight section at the nozzle exit
div_half_angle_2=divergent_half_angle # half angle of the diverging section at the termination
total_l=(L_cylindrical+Lconv+Ldiv) # total length of the engine
bell_r=Re # radius of the nozzle bell at the rim
domain_r=6*bell_r # radius of the area behind the nozzle for the plume
domain_ext=1.5*total_l+8 # extra length for the plume
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

# arc central point in the diverging bell - only if a bell arc is requested
if div_r is not None and div_r > 0:
    midangle = (div_half_angle_1 + div_half_angle_2) / 2.0
    x_center = points[14, 0] + div_r * math.sin(math.pi * div_half_angle_1 / 180.0)
    y_center = points[14, 1] - div_r * math.cos(math.pi * div_half_angle_1 / 180.0)
    aux_points[2] = (x_center - div_r * math.sin(math.pi * midangle / 180.0),
                     y_center + div_r * math.cos(math.pi * midangle / 180.0), 0)
else:
    # no bell arc for pure cone
    aux_points[2] = (math.nan, math.nan, 0)

#plt.scatter(aux_points[0:3,0], aux_points[0:3,1])
plt.show()

captured_output = io.StringIO()  # Create a StringIO object to capture the output


# convert points to wedge

def rotate_points(points,angle):
    for point in points:
        point[1]=point[1]*math.cos(math.pi*angle/180)
        point[2]=point[1]*math.sin(math.pi*angle/180)
    return points

# generate main and aux wedge points by rotation
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
    print(line)                     # keep console output
    vertices_lines.append(line)

print(");")

# now create the vertices_text used in replacements
vertices_text = "\n".join(vertices_lines)


axis_point_numbers=np.arange(0,9) #  poit numbers on the axis
wedge_point_numbers_1=np.arange(9,20) #  point numbers on the first wedge
wedge_point_numbers_2=np.arange(20,31) #  point numbers on the second wedge
# number of generated blocks
n_blocks=7

# define function that creates pseudo 4-point cross-sections for the main engine part
def contour(i):
    return [axis_point_numbers[i], wedge_point_numbers_2[i], wedge_point_numbers_1[i], axis_point_numbers[i]]

contours=[contour(i) for i in range(n_blocks+1)] #  define list of contours
blocks=[contours[i]+contours[i+1] for i in range(7)] #define list of blocks

y_cells=[20 for i in range(n_blocks)] # define list of cell numbers in y-direction  (uniform)
x_cells=[int((points[i+1,0]-points[i,0])/0.005+1) for i in range(n_blocks)] # define list of cell numbers in x-direction  (every 5 mm)

# generate block section of blockmeshdict
print("blocks")
print("(")

blocks_text = []
for i,block in enumerate(blocks):
    
    block = (f'\t hex ({block[0]} {block[1]} {block[2]} {block[3]} {block[4]} {block[5]} {block[6]} {block[7]} ) ({y_cells[i]} 1 {x_cells[i]}) simpleGrading (1 1 1)')
    print(block)
    blocks_text.append(block)

    blocks_text2 = "\n".join(blocks_text)
        


# number of generated blocks

# generate asym1 surfaces for wedge BC
asym1=[[contours[i][0],contours[i][2],contours[i+1][2],contours[i+1][0]] for i in range(n_blocks) ]
print (f'asym1 patches')

asym_text22 = []
for patch in asym1:
    asym1_text = (f'\t ({patch[0]} {patch[1]} {patch[2]} {patch[3]})')
    print(asym1_text)
    asym_text22.append(asym1_text)
    asym1_text2 = "\n".join(asym_text22)

# generate asym2 surfaces for wedge BC

asym2=[[contours[i+1][0],contours[i+1][1],contours[i][1],contours[i][0]] for i in range(n_blocks) ]
print (f'asym2 patches')
asym_text = []
for patch in asym2:
    asym_text2 = (f'\t({patch[0]} {patch[1]} {patch[2]} {patch[3]})')
    print(asym_text2)
    asym_text.append(asym_text2)
    asym2_text = "\n".join(asym_text)
    
# generate nozzle surfaces for nozzle wall BC
#nozzle[i]=np.array((wedges[i,2],wedges[i,1], wedges[i+1,1], wedges[i+1,2]))
nozzle_text = []
nozzle=[[contours[i][2],contours[i][1],contours[i+1][1],contours[i+1][2]] for i in range(n_blocks) ]
print (f'nozzle patches')
for patch in nozzle:
    nozzle2 = (f'\t({patch[0]} {patch[1]} {patch[2]} {patch[3]})')
    print(nozzle2)
    nozzle_text.append(nozzle2)
    nozzle_text2 = "\n".join(nozzle_text)

# generate inlet surfaces for inlet BC
print (f'inlet patch')
inlet=contours[0]
inlet_text = []
inlet2 = (f'\t({inlet[0]} {inlet[1]} {inlet[2]} {inlet[3]} )')
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

# generate outlet surfaces for outlet BC

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
print(edge2_2)
edge2_text2 = "\n".join(edge2_text)




if div_r is not None and div_r > 0:
    """
    Prints the arc segment for the nozzle extension bell if div_r is specified by the user based on their bell geometry
    """
    print(f'\t arc 14 15 ({aux_points[2,0]} {aux_points[2,1]} {aux_points[2,2]})')
    #edge5_text = captured_output.getvalue()


"""
Prints the arc segments at the throat section of the nozzle
"""
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
for old, new in replacements.items(): #this replaces the ambient pressure and chamber pressure specified by the user into the pressure file for OpenFOAM
    blockMeshDict_text = blockMeshDict_text.replace(old, new)

block_file = open(blockMeshDict, 'w')
block_file.write(blockMeshDict_text)
block_file.close()

os.makedirs(folder_name, exist_ok=True)
os.makedirs(f"{folder_name}/0", exist_ok=True)
os.makedirs(f"{folder_name}/constant", exist_ok=True)
os.makedirs(f"{folder_name}/system", exist_ok=True)


for filename in ["blockMeshDict.txt"]:
    new_name = os.path.splitext(filename)[0]
    shutil.move(filename, os.path.join(folder_name, "system", new_name))

os.chdir(folder_name)
