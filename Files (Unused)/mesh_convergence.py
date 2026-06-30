import sys
from matplotlib import pyplot as plt
import numpy as np

refinement_ratios = [1,4,16.015,32.0302 ]

#mach numbers at 500 in excel
#[3.6732898,3.6282701,3.6389852]
#Exit Mach ($M_e$)




independent_var = [2812.60,2879.50,2927.86, 2936.872944] #axial velocity at nozzle exit



plt.figure()
plt.plot(refinement_ratios, independent_var, 'o-', color='black', markersize=2)

plt.xlabel('Refinement Ratio [$N_{coarse}$/$N_{mesh}$]'   , fontsize=14)
plt.ylabel('$U_x$ [m/s]', fontsize=14)
plt.grid()
plt.show()