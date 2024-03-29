# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 15:14:43 2022

@author: benja
"""

import numpy as np
import matplotlib.pyplot as plt

# functions ----------------------------------------
def temp_to_voltage(temp_in, T_ref): 
    # convert temperature (C) to EMF (mV) type T thermocouple+
    # only for -270C <= temp_in <= 400C
    emf_pre = 0 # preliminary voltage before subtracting out emf_ref
    for ind in range(len(b[temp_in>0])):
        emf_pre += b[temp_in>0][ind]*(temp_in**ind)
    emf_ref = 0
    for ind in range(len(b[T_ref>0])):
        emf_ref += b[T_ref>0][ind]*(T_ref**ind)
    return emf_pre - emf_ref # returns emf in millivolts (mV)

def voltage_to_temp(emf_measured, T_ref):
    # convert EMF (mV) to temperature (C) for type T thermocouple
    # only for -5.603mV <= emf_measured <= 20.872mV
    temp_out = 0
    emf_ref = 0 # compute EMF for reference temperature using polynomial
    for ind in range(len(b[T_ref>0])):
        emf_ref += b[T_ref>0][ind]*(T_ref**ind)
    final_emf = emf_measured + emf_ref
    for ind in range(len(c[final_emf>0])):
        temp_out += c[final_emf>0][ind]*(final_emf**ind)
    return temp_out # returns temperature in Celsius (C)

def kelvin_to_celsius(kelv_in):
    return kelv_in - 273.15

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

def invoke_plot_params(filename="unnamed_1200dpi"):
    plt.tick_params(axis='both',which='major', direction = 'in', 
                    top = 1, right = 1, length=6,width=1,labelsize=16)
    plt.tick_params(axis='both',which='minor', direction = 'in', 
                    top = 1, right = 1, length=2,width=1,labelsize=16)
    plt.grid()
    # plt.savefig("plot_img_1200dpi/"+filename, dpi=1200, bbox_inches='tight')
    plt.show()

def plot_polynomials(temp_range, volt_range, Tref):
    # temp_range units: C   volt_range units: mV   Tref units: C
    temps_in = np.linspace(temp_range[0], temp_range[1], 301) # units: C
    volts_out = [temp_to_voltage(temp, Tref_C) for temp in temps_in]
    
    plt.plot(temps_in, volts_out)
    plt.title('Temperature to Voltage Conversion Polynomial', pad=20)
    plt.xlabel('Temperature (C)')
    plt.ylabel('Voltage (mV)')
    invoke_plot_params("temp_to_voltage")
    
    volts_in = np.linspace(volt_range[0], volt_range[1], 301) # units: mV
    temps_out = [voltage_to_temp(volt, Tref_C) for volt in volts_in]
    
    volt_diff = 0.05 # mV
    cold_start = 0.03153 # mV
    cold_volts = [cold_start, cold_start + volt_diff]
    hot_start = 0.06355 # mV
    hot_volts = [hot_start, hot_start + volt_diff]
    cold1 = voltage_to_temp(cold_volts[0], Tref_C)
    cold2 = voltage_to_temp(cold_volts[1], Tref_C)
    hot1 = voltage_to_temp(hot_volts[0], Tref_C)
    hot2 = voltage_to_temp(hot_volts[1], Tref_C)
    # volts = [0.00697,
    #         0.10697,
    #         0.06355,
    #         0.16355,
    #         0.00348,
    #         0.00348,
    #         0.03153,
    #         0.03153]
    # temps = []
    # for volt in volts:
    #     temps.append(voltage_to_temp(volt, Tref_C))
    
    plt.plot(volts_in, temps_out, 'k')
    plt.plot(cold_volts, [cold1, cold2], 'bo',
              label='Temp diff: %.2f C' % (cold2-cold1))
    plt.plot(hot_volts, [hot1, hot2], 'ro',
              label='Temp diff: %.2f C' % (hot2-hot1))
    # plt.plot(volts[0:4], temps[0:4], 'ro')
    # plt.plot(volts[0:4], temps[0:4], 'ro')
    # plt.plot(volts[4:], temps[4:], 'bo')
    # plt.plot(volts[4:], temps[4:], 'bo')
    # plt.axvline(0, -195, -180, color='r') # for plotting horizontal lines
    # test = plt.axvline(x=.2, ymin=-195, ymax=-180)
    # plt.plot(test)
    plt.title("Voltage to Temperature Conversion Polynomial, Tref=%d K" %
              (T_ref_K), pad=20)
    plt.xlabel('Voltage (mV)', fontsize=font_size)
    plt.ylabel('Temperature (C)', fontsize=font_size)
    plt.legend(bbox_to_anchor=(1.05,1))
    invoke_plot_params("voltage_to_temp")

font_size = 16
T_ref_K = 80 # units: K
Tref_C = kelvin_to_celsius(T_ref_K) # reference temperature in celsius
plot_polynomials([-200, 100], [-4, 20], Tref_C) # wide ranges
plot_polynomials([-200, -180], [-0.1, 0.3], Tref_C) # realistic ranges













