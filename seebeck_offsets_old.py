# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 22:03:37 2021

@author: benja
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 21:15:09 2021

@author: benja
"""
# OLD VERSION with writing data to text file
# functions ----------------------------------------
def temp_to_voltage(temp): 
    # convert temperature to EMF type T thermocouple
    emf_voltage = 0 
    for ind in range(len(b)):
        emf_voltage += b[ind]*(temp**ind)
    return emf_voltage

def voltage_to_temp(emf_voltage):
    # convert EMF to temperature for type T thermocouple
    temp = 0
    for ind in range(len(c)):
        temp += c[ind]*(emf_voltage**ind)
    return temp

def round_and_print(message, list_in, digits):
    print(message)
    list_to_print = [round(item, digits) for item in list_in]
    print(list_to_print)
    print('') # to get a new line

# main script ---------------------------------------

# import matplotlib.pyplot as plt
# dQ/dt = P = kA(dT/dx)
# power delivered to the sample
# powers = [0.5, 0.8, 1, 2, 5, 10, 20]  
powers = [0, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, \
          6e-3, 7e-3, 8e-3, 9e-3, 10e-3] # units: W NEED 100 values here
round_and_print("Input power values: ", powers, 7)
# thermal conductivity kappa for Bi2Te3
kappa = 2.0 # units: W/(m*K)
# 2mm x 2mm cross sectional area 
area = 0.002**2 # units: m^2
# location of hot and cold thermocoulpe probe points
x_cold = 0.001667 # 1.67mm, or sample length / 3
# length difference between each thermocouple probe point
dx = 0.001667 # 5mm
x_hot = x_cold + dx # 30mm
# difference in temperature between each thermocouple (delta T)
dT = [(pwr * dx)/(kappa * area) for pwr in powers]
# derivative of T with respect to x is dT/dx (slope)
# T(x) = (dT/dx)x + Tref
# reference temperature at the base of the sample
Tref = 80 # units: K

Thots = [(delta_T/dx)*x_hot + Tref for delta_T in dT]
Tcolds = [(delta_T/dx)*x_cold + Tref for delta_T in dT]
# get difference in temperature (dT and delta_temps are equal)
# delta_temps = [Thots[ind]-Tcolds[ind] for ind in range(len(Thots))]    

round_and_print("Initial temperature differences: ", dT, 7)

# introduce simulated temperature offsets here

# coefficients for EMF as a function of temperature b0-b14 (-270C<=T<=0C)
# NEED TO ACCOUNT FOR REFERENCE TEMP ON ALL CONVERSIONS
# POLYNOMIALS ASSUME 0K REF TEMP 
# FIGURE OUT HOW TO CONVERT TEMP TO VOLTAGE
# SEE SUMMARY ON PAGE 11 FOR CONVERTING VOLTAGE TO TEMP (CODE ALREADY WRITTEN)
# coefficients for EMF as a function of temperature b0-b14 (-270C<=T<=0C)
b = [0, 3.8748106364e-2, 4.4194434347e-5, 1.1844323105e-7, # temp -> voltage
     9.0138019559e-10, 2.2651156593e-11, 3.6071154205e-13,
     3.8493939883e-15, 3.8493939883e-15, 2.8213521925e-17,
     1.4251594779e-19, 4.8768662286e-22, 1.0795539270e-24,
     1.3945027062e-27, 7.9795153927e-31]
# coefficients for temp to voltage for positive temp (0C<=T<=400C)
b_pos = [0, 3.8748106364e-2, 3.3292227880e-5, 2.0618243404e-7,
         -2.1882256846e-9, 1.0996880928e-11, -3.0815758772e-14,
         4.5479135290e-17, -2.7512901673e-20]

# coefficients for temperature as a function of EMF c0-7 (-200C<=T<=0C)
c = [0, 2.5949192e1, -2.1316967e-1, 7.9018692e-1, # voltage -> temp 
     4.2527777e-1, 1.3304473e-1, 2.0241446e-2,
     1.2668171e-3]
# coefficients for voltage to temp for positive temp/votlage (0C<=T<=400C)
c_pos = [0, 2.592800e1, -7.602961e-1, 4.637791e-2, 
         2.165394e-3, 6.048144e-5, -7.293422e-7, 0]
    
# use conversion polynomials to get delta V values
delta_V34 = [temp_to_voltage(temp) for temp in Thots]
delta_V12 = [temp_to_voltage(temp) for temp in Tcolds]
round_and_print("Voltage across hot thermocouple: ", delta_V34, 7)
round_and_print("Voltage across cold thermocouple: ", delta_V12, 7)

# plt.plot(delta_V34)

# use polynomials to return 
new_Thots = [voltage_to_temp(volt) for volt in delta_V34]
new_Tcolds = [voltage_to_temp(volt) for volt in delta_V12]
round_and_print("New hot temperatures: ", new_Thots, 7)
round_and_print("New cold temperatures: ", new_Tcolds, 7)

data_file = open("data_file.txt", "a+")
data_file.truncate(0) # clear file contents
for item in delta_V34:
    print(item, file=data_file)
print('', file=data_file)
for item in delta_V12:
    print(item, file=data_file)
for item in b:
    print(item, file=data_file)
for item in c:
    print(item, file=data_file)
    
    
data_file.close()

new_dT = [new_Thots[ind]-new_Tcolds[ind] for ind in range(len(new_Thots))]

# NIST provided Seebeck coefficient for Bi2Te3+x
S_nist = -74.1 # units: uV/K

deltaV_seebecks = [-1*S_nist * delta_T for delta_T in new_dT]






















