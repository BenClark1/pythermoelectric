# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 16:59:03 2021

@author: benja
"""
import matplotlib.pyplot as plt
import math

def Seebeck_Cu(T):
# T must be in kelvin: Seebeck coeff [uV/K]
    exp_term = math.exp(-1*T/93)
    rational_fraction = 0.442/(1+(T/172.4)**3)
    s_cu = 0.041*T*( exp_term + 0.123 -  rational_fraction) + 0.804
    return s_cu

def Seebeck_SRM3451(T):
# Standard Reference Material (SRM) 3451 = Bi2Te3+x
# NIST provided Seebeck coefficients for SRM3451
# T must be in kelvin: Seebeck coeff [uV/K]
    A = 295 # entral temp (room temp) K    
    S_A = -230.03 # µV/K
    a = -2.2040e-1 # coefficients
    b = 3.9706e-3
    c = 7.2922e-6 
    d = -1.0864e-9     
    term1 = a*T*(1-(A/T))
    term2 = b*T**2*(1-(A/T))**2
    term3 = c*T**3*(1-(A/T))**3
    term4 = d*T**4*(1-(A/T))**4
    return S_A + term1 + term2 + term3 + term4
    
def Seebeck_constantan(T): # units of T: K
    # get the Seebeck coefficent of Constantan
    # first get Seebeck coefficient for Copper/Constantan T-type thermocouple
    # MUST BE FOR 73.15 K < T < 673.15 K
    S_Cu_Con = 4.37184 + 0.1676*T - (1.84371e-4)*T**2 \
        + (1.2244e-7)*T**3 - (4.47618e-11)*T**4
    # subtract Cu/Con Seebeck value from copper Seebeck to get Con Seebeck
    return Seebeck_Cu(T) - S_Cu_Con

powers = [0, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, \
          6e-3, 7e-3, 8e-3, 9e-3, 10e-3] # units: W
kappa = 2.0 # units: W/(m*K)  thermal conductivity kappa for Bi2Te3
area = 0.002**2 # units: m^2  2mm x 2mm cross sectional area 
x_cold = 0.001667 # 1.67mm, or sample length / 3
dx = 0.001667 # 1.67mm  length difference between each thermocouple probe point
x_hot = x_cold + dx # 2*1.67mm
# difference in temperature between each thermocouple (delta T)
dT_true = [(pwr * dx)/(kappa * area) for pwr in powers] #units: K
T_ref_K = 80 # units: K
Thots = [(delta_T/dx)*x_hot + T_ref_K for delta_T in dT_true]
Tcolds = [(delta_T/dx)*x_cold + T_ref_K for delta_T in dT_true] 

avg_temps = [(Thots[i]+Tcolds[i])/2 for i in range(len(Thots))]
changing_seebecks = [Seebeck_SRM3451(avg_temp) for avg_temp in avg_temps]
control = [Seebeck_SRM3451(T_ref_K)] * len(changing_seebecks)

plt.plot(avg_temps, changing_seebecks, label='Using avg. temp.')
plt.plot(avg_temps, control, label='Using reference temp.')
plt.title('Seebeck coefficient calculated with reference temperature and average temperature', pad=20)
plt.xlabel(r'Average Temperature (K)')
plt.ylabel('Seebeck Coefficient (uV/K)')
plt.legend(bbox_to_anchor=(1.05,1))
# plt.autoscale(enable=False, axis='y')
plt.grid()
plt.show()