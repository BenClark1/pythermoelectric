# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 21:15:09 2021

@author: benja
"""
import matplotlib.pyplot as plt
import math
# functions ----------------------------------------
def temp_to_voltage(temp_in, T_ref): 
    # convert temperature (C) to EMF (mV) type T thermocouple
    emf_pre = 0 # preliminary voltage before subtracting out emf_ref
    for ind in range(len(b[temp_in>0])):
        emf_pre += b[temp_in>0][ind]*(temp_in**ind)
    emf_ref = 0
    for ind in range(len(b[T_ref>0])):
        emf_ref += b[T_ref>0][ind]*(T_ref**ind)
    return emf_pre - emf_ref # returns emf in millivolts (mV)

def voltage_to_temp(emf_measured, T_ref):
    # convert EMF (mV) to temperature (C) for type T thermocouple
    temp_out = 0
    emf_ref = 0 # compute EMF for reference temperature using polynomial
    for ind in range(len(b[T_ref>0])):
        emf_ref += b[T_ref>0][ind]*(T_ref**ind)
    final_emf = emf_measured + emf_ref
    for ind in range(len(c[final_emf>0])):
        temp_out += c[final_emf>0][ind]*(final_emf**ind)
    return temp_out # returns temperature in Celsius (C)

def round_and_print(message, list_in, digits):
    print(message)
    list_to_print = [round(item, digits) for item in list_in]
    print(list_to_print)
    print('') # to get a new line
    
def kelvin_to_celsius(kelv_in):
    return kelv_in - 273.15

def calculate_trendline(x_vals, y_vals):
    # for ind in range(len(x_vals)):
    #     # current_ind = x_vals.index(uvolt)
    #     print(x_vals[ind], end='')
    #     print('\t', end='')
    #     print(y_vals[ind])
    
    mean_x = sum(x_vals) / len(x_vals)
    mean_y = sum(y_vals) / len(y_vals)
    
    variations_x = [(xi-mean_x)**2 for xi in x_vals]
    variations_y = [(yi-mean_y)**2 for yi in y_vals]
    
    # Standard Deviations for x and y
    sdev_x = math.sqrt(sum(variations_x)/len(x_vals))
    sdev_y = math.sqrt(sum(variations_y)/len(y_vals))
    
    correlation = 0
    for ind in range(len(x_vals)):
        correlation += (x_vals[ind] - mean_x)*(y_vals[ind] - mean_y)
    correlation /= math.sqrt(sum(variations_x))
    correlation /= math.sqrt(sum(variations_y))
    
    trend_slope = correlation * (sdev_y/sdev_x)
    trend_intercept = mean_y - (trend_slope * mean_x)
    trend_y_vals = [trend_slope*xval + trend_intercept for xval in x_vals]
    
    # print("------")
    # print("Calculated trendline: ")
    # for ind in range(len(x_vals)):
    #     print(x_vals[ind], end='')
    #     print('\t', end='')
    #     print(trend_y_vals[ind])
        
    print("\nTrendline slope: %f" % (trend_slope))
    print("Trendline intercept: %f" % (trend_intercept))   
    print("Correlation coefficient: %f" % (correlation)) 

    return {'slope':trend_slope, 
            'intercept':trend_intercept, 
            'trendline':trend_y_vals}

def Seebeck_Cu(Tref_K):
# Tref must be in kelvin, gives Roberts data: Seebeck coeff [uV/K]
    exp_term = math.exp(-1*Tref_K/93)
    rational_fraction = 0.442/(1+(Tref_K/172.4)**3)
    s_cu = 0.041*Tref_K*( exp_term + 0.123 -  rational_fraction) + 0.804
    return s_cu

def get_s_coeff(Tref_K):
    # NIST provided Seebeck coefficients for Bi2Te3+x
    # keys are ref temps in K, values are Seebeck coefficients
    S_nist = {30.11: -33.36, # units: uV/K
            40.14: -40.42,
            50.16: -48.01,
            60.25: -56.26,
            70.29: -65.28,
            80.33: -74.1,
            100.34: -92.79,
            120.37: -111.6,
            140.38: -129.6,
            160.4: -147.22,
            180.41: -163.81,
            200.43: -179.4,
            220.49: -193.96,
            240.52: -206.22,
            260.52: -217.26,
            280.71: -226.11,
            300.73: -231.36,
            310.74: -232.82 }
    diffs = [abs(list(S_nist.keys())[ind]-Tref_K) for ind in range(len(S_nist))]
    coeff = None
    for key in S_nist.keys():
        if (round(abs(key-Tref_K),5) == round(min(diffs),5)):
            coeff = S_nist[key]
    return coeff
    

# main script ---------------------------------------

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
dx = 0.001667 # 1.67mm
x_hot = x_cold + dx # 2*1.67mm
# difference in temperature between each thermocouple (delta T)
dT = [(pwr * dx)/(kappa * area) for pwr in powers] #units: K
# derivative of T with respect to x is dT/dx (slope)
# T(x) = (dT/dx)x + T_ref_K
# reference temperature at the base of the sample
# T_ref_K = 273.15 # units: K
T_ref_K = 80 # units: K
# get hot and cold temperatures and convert to celsius
Thots_C = [kelvin_to_celsius((delta_T/dx)*x_hot + T_ref_K) for delta_T in dT]
Tcolds_C = [kelvin_to_celsius((delta_T/dx)*x_cold + T_ref_K) for delta_T in dT] 
T_ref_C = kelvin_to_celsius(T_ref_K) # reference temperature in celsius

# round_and_print("Initial temperature differences: ", dT, 7)
round_and_print("Initial hot temperatures (Thot C): ", Thots_C, 7)
round_and_print("Initial cold temperatures(Tcold C): ", Tcolds_C, 7)

# coefficients for EMF as a function of temperature b0-b14 (-270C<=T<=0C)
b_neg = [0.0,
0.387481063640e-1,
0.441944343470e-4,
0.118443231050e-6,
0.200329735540e-7,
0.901380195590e-9,
0.226511565930e-10,
0.360711542050e-12,
0.384939398830e-14,
0.282135219250e-16,
0.142515947790e-18,
0.487686622860e-21,
0.107955392700e-23,
0.139450270620e-26,
0.797951539270e-30]
# coefficients for temp to voltage for positive temp (0C<=T<=400C)
b_pos = [0.0,
0.387481063640e-1,
0.332922278800e-4,
0.206182434040e-6,
-0.218822568460e-8,
0.109968809280e-10,
-0.308157587720e-13,
0.454791352900e-16,
-0.275129016730e-19]
b = [b_neg, b_pos] # b[1] gives positive polynomial, b[0] gives negative
# coefficients for temperature as a function of EMF c0-7 (-200C<=T<=0C)
c_neg = [0.0,
2.5949192e1,
-2.1316967e-1,
7.9018692e-1,
4.2527777e-1,
1.3304473e-1,
2.0241446e-2,
1.2668171e-3] 
# coefficients for voltage to temp for positive temp/voltage (0C<=T<=400C)
c_pos = [0.0,
2.592800e1,
-7.602961e-1,
4.637791e-2,
-2.165394e-3,
6.048144e-5,
-7.293422e-7]
c = [c_neg, c_pos] # c[1] gives positive polynomial, c[0] gives negative
# use conversion polynomials to get delta V values
delta_V34 = [temp_to_voltage(temp, T_ref_C) for temp in Thots_C]
delta_V12 = [temp_to_voltage(temp, T_ref_C) for temp in Tcolds_C]
round_and_print("Voltage across hot thermocouple: ", delta_V34, 7)
round_and_print("Voltage across cold thermocouple: ", delta_V12, 7)

# introduce simulated voltage offsets here

# use polynomials to return to temperatures
new_Thots_C = [voltage_to_temp(volt, T_ref_C) for volt in delta_V34]
new_Tcolds_C = [voltage_to_temp(volt, T_ref_C) for volt in delta_V12]
round_and_print("New hot temperatures (C): ", new_Thots_C, 7)
round_and_print("New cold temperatures(C): ", new_Tcolds_C, 7)
# new_dT is the same in both Kelvin and Celsius
new_dT = [new_Thots_C[ind]-new_Tcolds_C[ind] for ind in range(len(new_Thots_C))]

S_Cu = round(Seebeck_Cu(T_ref_K), 3) # units: uV/K
# deltaV_seebecks = [-1*S_nist * delta_T for delta_T in new_dT] #not caclulated

true_deltaV = [-1*(get_s_coeff(T_ref_K) - S_Cu)*delta_T for delta_T in new_dT]
# note: true_deltaV is in uV
# get a dictionary with slope, intercept, and trendline y values
trend_info = calculate_trendline(new_dT, true_deltaV)


plt.plot(new_dT, true_deltaV, 'r.', new_dT, trend_info['trendline'], 'b')
plt.title('Thermoelectric Votlage Produced by Seebeck Effect in Bi₂Te₃₊ₓ')
plt.xlabel('Temperature Difference (K)')
plt.ylabel('Thermoelectric Voltage (uV)')
plt.show()

S_sample = -1*trend_info['slope'] + S_Cu
print("\nFinal Seebeck Coefficient of the Sample: ")
print(round(S_sample, 7))

print("\nDifferences between original and new temps: ")
print("hot: ")
print([Thots_C[i]-new_Thots_C[i] for i in range(len(new_Thots_C))])
print("cold: ")
print([Tcolds_C[i]-new_Tcolds_C[i] for i in range(len(new_Tcolds_C))])







