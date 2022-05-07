# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 11:50:54 2022

@author: benja
"""

import numpy as np
import matplotlib.pyplot as plt
# from seebeck_offset_sim import (temp_to_voltage, voltage_to_temp, 
#                                 kelvin_to_celsius, calculate_trendline, 
#                                 Seebeck_platinum, invoke_plot_params)
from seebeck_offset_sim import (Seebeck_platinum, invoke_plot_params)

temp_range = [100, 600]
num_points = 51
T = np.linspace(temp_range[0], temp_range[1], num_points)
S_Pt = [Seebeck_platinum(temp) for temp in T]

plt.plot(T, S_Pt)
plt.title('Temperature Dependent Seebeck Coefficient of Platinum', pad=20)
plt.xlabel('Temperature (K)')
plt.ylabel('Seebeck Coefficient of Platinum (uV / K)')
plt.xlim([50, 650])
plt.ylim([-15, 10])
invoke_plot_params("S_Pt_vs_temp")
