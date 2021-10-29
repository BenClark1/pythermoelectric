# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 21:15:09 2021

@author: benja
"""
import matplotlib.pyplot as plt
import get_slope as gs
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
# T(x) = (dT/dx)x + T_reference
# reference temperature at the base of the sample
# T_reference = 273.15 # units: K
T_reference = 80 # units: K
# get hot and cold temperatures and convert to celsius
Thots_C = [kelvin_to_celsius((delta_T/dx)*x_hot + T_reference) for delta_T in dT]
Tcolds_C = [kelvin_to_celsius((delta_T/dx)*x_cold + T_reference) for delta_T in dT] 
T_ref_C = kelvin_to_celsius(T_reference) # reference temperature in celsius

# round_and_print("Initial temperature differences: ", dT, 7)
round_and_print("Initial hot temperatures (Thot C): ", Thots_C, 7)
round_and_print("Initial cold temperatures(Tcold C): ", Tcolds_C, 7)

# OLD VALUES
# # coefficients for EMF as a function of temperature b0-b14 (-270C<=T<=0C)
# b_neg = [0, 3.8748106364e-2, 4.4194434347e-5, 1.1844323105e-7, # temp -> voltage
#       9.0138019559e-10, 2.2651156593e-11, 3.6071154205e-13,
#       3.8493939883e-15, 3.8493939883e-15, 2.8213521925e-17,
#       1.4251594779e-19, 4.8768662286e-22, 1.0795539270e-24,
#       1.3945027062e-27, 7.9795153927e-31]
# # coefficients for temp to voltage for positive temp (0C<=T<=400C)
# b_pos = [0, 3.8748106364e-2, 3.3292227880e-5, 2.0618243404e-7,
#           -2.1882256846e-9, 1.0996880928e-11, -3.0815758772e-14,
#           4.5479135290e-17, -2.7512901673e-20]

# # coefficients for temperature as a function of EMF c0-7 (-200C<=T<=0C)
# c_neg = [0, 2.5949192e1, -2.1316967e-1, 7.9018692e-1, # voltage -> temp 
#       4.2527777e-1, 1.3304473e-1, 2.0241446e-2,
#       1.2668171e-3]
# # coefficients for voltage to temp for positive temp/votlage (0C<=T<=400C)
# c_pos = [0, 2.592800e1, -7.602961e-1, 4.637791e-2, 
#           2.165394e-3, 6.048144e-5, -7.293422e-7, 0]
    
# NEW VALUES
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
# coefficients for voltage to temp for positive temp/votlage (0C<=T<=400C)
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

# plt.plot(delta_V34)
# introduce simulated voltage offsets here
# use polynomials to return to temperatures

new_Thots_C = [voltage_to_temp(volt, T_ref_C) for volt in delta_V34]
new_Tcolds_C = [voltage_to_temp(volt, T_ref_C) for volt in delta_V12]
round_and_print("New hot temperatures (C): ", new_Thots_C, 7)
round_and_print("New cold temperatures(C): ", new_Tcolds_C, 7)

new_dT = [new_Thots_C[ind]-new_Tcolds_C[ind] for ind in range(len(new_Thots_C))]

# NIST provided Seebeck coefficient for Bi2Te3+x
S_nist = -74.1 # units: uV/K
S_Cu = 6.5 # units: uV/K
deltaV_seebecks = [-1*S_nist * delta_T for delta_T in new_dT]

true_deltaV = [-1*(S_nist - S_Cu)*delta_T for delta_T in new_dT]
# get a dictionary with slope, intercept, and trendline y values
trend_info = gs.calculate_trendline(new_dT, true_deltaV)


plt.plot(new_dT, true_deltaV, 'ro', new_dT, trend_info['trendline'], 'b')
plt.xlabel('new_dT')
plt.ylabel('true_deltaV')
plt.show()

S_sample = -1*trend_info['slope'] + S_Cu
print("\nFinal Seebeck Coefficient of the Sample: ")
print(round(S_sample, 7))

# plotting stuff
# import matplotlib.pyplot as plt
# plt.semilogx(ConcVals1,SeebeckVals1/10,'--b')
# plt.semilogx(ConcVals2,SeebeckVals2/10,'--r')
# plt.semilogx(ConcVals3,SeebeckVals3/10,'--k')
# plt.semilogx(ConcValsnH,SeebeckVals3/10,'--m')
# plt.title('Witting Figure 10.A: Varying methods')
# plt.xlabel('n [cm$^{-3}$]')
# plt.ylabel('|S| [$\mu$V/K]')
# plt.legend(['Generalized Fermi','Generalized Transport','j$^{th}$ order Fermi','j$^{th}$ order Fermi and n$_{H}$'], loc = 'upper right')
# plt.show()


# test cases tell whether or not these values changed from previous version
# new_colds_temp = [-193.1568675, -993.3270197, -992.8297834, 
#                   -992.3319282, -991.8334559, -991.3343678, -990.8346655, 
#                   -990.3343504, -989.8334241, -989.331888, -988.8297436]
# delta_V12_temp = [0.0, 0.0034801, 0.0069661, 0.0104581, 0.0139559, 
#                   0.0174597, 0.0209694, 0.0244851, 0.0280066, 0.0315341, 
#                   0.0350674]
# print(new_colds_temp==[round(item, 7) for item in new_Tcolds_C])
# print(delta_V12_temp==[round(item, 7) for item in delta_V12])

print("\nError values between original and new temps: ")
print("hot: ")
print([Thots_C[i]-new_Thots_C[i] for i in range(len(new_Thots_C))])
print("cold: ")
print([Tcolds_C[i]-new_Tcolds_C[i] for i in range(len(new_Tcolds_C))])







