# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 21:15:09 2021

@author: benja
"""
import matplotlib.pyplot as plt
import math
import numpy as np

# functions ----------------------------------------
def temp_to_voltage(temp_in, T_ref): 
    # convert temperature (C) to EMF (mV) type T thermocouple
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

def round_and_print(message, list_in, digits):
    print(message)
    list_to_print = [round(item, digits) for item in list_in]
    print(list_to_print)
    print('')
    
def kelvin_to_celsius(kelv_in):
    return kelv_in - 273.15

def calculate_trendline(x_vals, y_vals):

    if len(x_vals) != len(y_vals):
        print("Trendline failed: lists must have equal length")
        return None
    
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
        
    # print("\nTrendline slope: %f" % (trend_slope))
    # print("Trendline intercept: %f" % (trend_intercept))   
    # print("Correlation coefficient: %f" % (correlation)) 

    return {'slope':trend_slope, 
            'intercept':trend_intercept, 
            'trendline':trend_y_vals}

def Seebeck_Cu(T):
# T must be in kelvin: Seebeck coeff [uV/K]
    exp_term = math.exp(-1*T/93)
    rational_fraction = 0.442/(1+(T/172.4)**3)
    s_cu = 0.041*T*( exp_term + 0.123 -  rational_fraction) + 0.804
    return s_cu

def Seebeck_SRM3451(T):
# Standard Reference Material (SRM) 3451 = Bi2Te3+x
# NIST provided Seebeck coefficients for SRM3451
# T must be in KELVIN: Seebeck coeff [uV/K]
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

def seebeck_measurement(Thots_C, Tcolds_C, offs, plot=False, print_vals=False):
    # offs: a list of 5 offsets, first term must be zero
        # index of offs corresponds to offset location 
    # use conversion polynomials to get delta V values
    delta_V12_true = [temp_to_voltage(temp, Tref_C) for temp in Thots_C] # mV
    delta_V34_true = [temp_to_voltage(temp, Tref_C) for temp in Tcolds_C] # mV
    # add simulated voltage offsets, convert to mV
    delta_V12_meas = [volt + offs[1]/1000 + offs[2]/1000 for volt in delta_V12_true]
    delta_V34_meas = [volt + offs[3]/1000 + offs[4]/1000 for volt in delta_V34_true]
    # round_and_print("Voltage across hot thermocouple: ", delta_V12_true, 7)
    # round_and_print("Voltage across cold thermocouple: ", delta_V34_true, 7)
    if print_vals: # print two thermocouple voltages to see realistic example
        print("true hot voltage LP: ", round(delta_V12_true[1], 5))
        print("offset hot voltage LP: ", round(delta_V12_meas[1], 5))
        print("true cold voltage LP: ", round(delta_V34_true[1], 5))
        print("offset cold voltage LP: ", round(delta_V34_meas[1], 5))    
        print("true hot voltage: ", round(delta_V12_true[9], 5))
        print("offset hot voltage: ", round(delta_V12_meas[9], 5))
        print("true cold voltage: ", round(delta_V34_true[9], 5))
        print("offset cold voltage: ", round(delta_V34_meas[9], 5))
    # use polynomials to return to temperatures
    offs_Thots_C = [voltage_to_temp(volt, Tref_C) for volt in delta_V12_meas]
    offs_Tcolds_C = [voltage_to_temp(volt, Tref_C) for volt in delta_V34_meas]
    # round_and_print("Offset hot temperatures (C): ", offs_Thots_C, 7)
    # round_and_print("Offset cold temperatures(C): ", offs_Tcolds_C, 7)
    # meas_dT is the same in both Kelvin and Celsius
    meas_dT = [offs_Thots_C[ind]-offs_Tcolds_C[ind] for
              ind in range(len(offs_Thots_C))]
        
    S_Cu = round(Seebeck_Cu(T_ref_K), 3) # units: uV/K
    S_Con = round(Seebeck_constantan(T_ref_K), 3) # units: uV/K
    true_deltaV13 = [-1*(Seebeck_SRM3451(T_ref_K) - S_Con)*delta_T for
                     delta_T in dT_true]
    true_deltaV24 = [-1*(Seebeck_SRM3451(T_ref_K) - S_Cu)*delta_T for
                     delta_T in dT_true] #dT_true is a global variable
    # note: true_deltaV is in uV
    # introduce voltage offset for true_deltaV lists
    meas_deltaV13 = [volt + offs[1] + offs[3] for volt in true_deltaV13] # uV
    meas_deltaV24 = [volt + offs[2] + offs[4] for volt in true_deltaV24] # uV
    meas_deltaV = [meas_deltaV24, meas_deltaV13] # uV
    use_top_13_wires = False # choose between deltaV13 or deltaV24 for Seebeck voltage
    
    # get a dictionary with slope, intercept, and trendline y values
    trend_info = calculate_trendline(meas_dT, meas_deltaV[use_top_13_wires])
    
    if plot:
        plt.plot(meas_dT, meas_deltaV[use_top_13_wires], 'b.')
        plt.plot(meas_dT, trend_info['trendline'], 'b')
        plt.title('Thermoelectric Votlage Produced by Seebeck Effect in Bi₂Te₃₊ₓ', 
                  pad=20)
        plt.xlabel('Measured $\Delta$ Temperature (K)', fontsize=font_size)
        plt.ylabel('Measured Seebeck Voltage (uV)', fontsize=font_size)
        invoke_plot_params("Seebeck_Effect_for_Slope")
        
    S_sample = -1*trend_info['slope'] + [S_Cu, S_Con][use_top_13_wires] # need to add S_const
    # print("\nFinal Seebeck Coefficient of the Sample: ")
    # print(round(S_sample, 9))
    
    # print("\nDifferences between original and new temps (C): ")
    # print("hot: ")
    # print([Thots_C[i]-offs_Thots_C[i] for i in range(len(offs_Thots_C))])
    # print("cold: ")
    # print([Tcolds_C[i]-offs_Tcolds_C[i] for i in range(len(offs_Tcolds_C))])

    return S_sample

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
                    top = 1, right = 1, length=6,width=1,labelsize=font_size)
    plt.tick_params(axis='both',which='minor', direction = 'in', 
                    top = 1, right = 1, length=2,width=1,labelsize=font_size)
    plt.grid()
    # plt.savefig("plot_img_1200dpi/"+filename, dpi=1200, bbox_inches='tight')
    plt.show()
    
    
def plot_DT_offset(offs, Tref_K): # DT_offset can be plotted vs. DT_true
    Tref_C_loc = kelvin_to_celsius(Tref_K)
    # use conversion polynomials to get delta V values
    delta_V12_true = [temp_to_voltage(temp, Tref_C_loc) for temp in Thots_C] # mV
    delta_V34_true = [temp_to_voltage(temp, Tref_C_loc) for temp in Tcolds_C] # mV
    # add simulated voltage offsets, convert to mV
    delta_V12_meas = [volt + offs_inputs[1]/1000 + offs_inputs[2]/1000 
                      for volt in delta_V12_true] # mV
    delta_V34_meas = [volt + offs_inputs[3]/1000 + offs_inputs[4]/1000 
                      for volt in delta_V34_true] # mV
    
    # use polynomials to return to temperatures
    offs_Thots_C = [voltage_to_temp(volt, Tref_C_loc) for volt in delta_V12_meas]
    offs_Tcolds_C = [voltage_to_temp(volt, Tref_C_loc) for volt in delta_V34_meas]
    meas_dT = [offs_Thots_C[ind]-offs_Tcolds_C[ind] for
          ind in range(len(offs_Thots_C))]
    
    offs_Delta_T = [meas_dT[ind]-dT_true[ind] for ind in range(len(dT_true))]
    return offs_Delta_T

def plot_polynomials(temp_range, volt_range, Tref):
    # temp_range units: C   volt_range units: mV   Tref units: C
    temps_in = np.linspace(temp_range[0], temp_range[1], num_points) # units: C
    volts_out = [temp_to_voltage(temp, Tref_C) for temp in temps_in]
    
    plt.plot(temps_in, volts_out)
    plt.title('Temperature to Voltage Conversion Polynomial', pad=20)
    plt.xlabel('Temperature (C)')
    plt.ylabel('Voltage (mV)')
    invoke_plot_params("temp_to_voltage")
    
    volts_in = np.linspace(volt_range[0], volt_range[1], num_points) # units: mV
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
    
    plt.plot(volts_in, temps_out, 'k')
    plt.plot(cold_volts, [cold1, cold2], 'bo',
              label='Temp diff: %.3f C' % (cold2-cold1))
    plt.plot(hot_volts, [hot1, hot2], 'ro',
              label='Temp diff: %.3f C' % (hot2-hot1))
    # plt.axvline(0, -195, -180, color='r') # for plotting horizontal lines
    # test = plt.axvline(x=.2, ymin=-195, ymax=-180)
    # plt.plot(test)
    plt.title("Voltage to Temperature Conversion Polynomial, Tref=%d K" %
              (T_ref_K), pad=20)
    plt.xlabel('Voltage (mV)', fontsize=font_size)
    plt.ylabel('Temperature (C)', fontsize=font_size)
    plt.legend(bbox_to_anchor=(1.05,1))
    invoke_plot_params("voltage_to_temp")

# main script ---------------------------------------

# dQ/dT_true = P = kA(dT_true/dx) = power delivered to the sample
powers = [0, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, \
          6e-3, 7e-3, 8e-3, 9e-3, 10e-3] # units: W eventually test 100 values here
round_and_print("Input power values (W): ", powers, 7)
kappa = 2.0 # units: W/(m*K)  thermal conductivity kappa for Bi2Te3
area = 0.002**2 # units: m^2  2mm x 2mm cross sectional area 
# location of hot and cold thermocoulpe probe points
x_cold = 0.001667 # 1.67mm, or sample length / 3
dx = 0.001667 # 1.67mm  length difference between each thermocouple probe point
x_hot = x_cold + dx # 2*1.67mm
# difference in temperature between each thermocouple (delta T)
dT_true = [(pwr * dx)/(kappa * area) for pwr in powers] #units: K
# derivative of T with respect to x is 
    # dT_true/dx (slope)  T(x) = (dT_true/dx)x + T_ref_K
# reference temperature at the base of the sample
T_ref_K = 80 # units: K
# T_ref_K = 293 # units: K
# get hot and cold temperatures and convert to celsius
Thots = [(delta_T/dx)*x_hot + T_ref_K for delta_T in dT_true]
Tcolds = [(delta_T/dx)*x_cold + T_ref_K for delta_T in dT_true] 
Thots_C = [kelvin_to_celsius(temp) for temp in Thots]
Tcolds_C = [kelvin_to_celsius(temp) for temp in Tcolds] 
Tref_C = kelvin_to_celsius(T_ref_K) # reference temperature in celsius
# create offsets in uV
# offset_list1 = [-200,-100,-50,-20,-10,-5,-2,-1,-0.5,-0.2,-0.1,-0.05,-0.02,-0.01, 
#          0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200]
main_offset_list = [-200, -100, -50, 0, 50, 100, 200]
offset_list1 = main_offset_list
offset_list2 = main_offset_list 
offset_list3 = main_offset_list
offset_list4 = main_offset_list
ind_zero = len(offset_list1)//2 # get index of zero offset (center of list) 

font_size = 16 

offs_inputs = [0, 0, 0, 0, 0]
num_points = 301
# for plotting true value
true_seebeck = [Seebeck_SRM3451(T_ref_K)] * num_points
percent_error = 5
minus_percent = [S*(1-(percent_error/100)) for S in true_seebeck]
plus_percent = [S*(1+(percent_error/100)) for S in true_seebeck]

# hold offsets 3 and 4 constant while varying 1 and 2 ------------ FIGURE 1
offset_list1 = np.linspace(main_offset_list[0],main_offset_list[-1],num_points)
for ind in range(len(offset_list2)):
    s_coeffs = []
    offs_inputs[2] = offset_list2[ind]
    for offset1 in offset_list1:
        offs_inputs[1] = offset1
        s_coeffs.append(seebeck_measurement(Thots_C, Tcolds_C, offs_inputs))
    
    plt.plot(offset_list1, s_coeffs, 
             label=r'$\delta V2=%.2f uV$' % (round(offs_inputs[2], 2)))
    

plt.plot(offset_list1, true_seebeck, 'r--', label="True Seebeck Coefficient")
# plt.plot(offset_list1, minus_percent, 'k:', 
#          label='%d'%(percent_error)+"% Error Bounds")
# plt.plot(offset_list1, plus_percent, 'k:')
plt.fill_between(offset_list1, minus_percent, plus_percent, color='gray', 
                 alpha=.3, label="$\u00B1$%d"%(percent_error)+
                 "% Typical Experimental\nUncertainty")
plt.title("Offsets in Hot Thermocouple Only, Tref=%d K"%(T_ref_K), pad=20)
plt.xlabel('$\delta$V1 ($\mu$V)', fontsize=font_size)
plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=font_size)
plt.legend(bbox_to_anchor=(1.05,1))
# plt.autoscale(enable=False, axis='y')
invoke_plot_params("S_vs_dV1_dV2")

# hold offsets 1 and 2 constant while varying 3 and 4------------ FIGURE 2
offset_list3 = offset_list1 # offset_list3 becomes a linspace for plot
offset_list1 = main_offset_list # reset offset_list1
offs_inputs = [0, 0, 0, 0, 0]
for ind in range(len(offset_list4)):
    s_coeffs = []
    offs_inputs[4] = offset_list4[ind]
    for offset3 in offset_list3:
        offs_inputs[3] = offset3
        s_coeffs.append(seebeck_measurement(Thots_C, Tcolds_C, offs_inputs))
    
    plt.plot(offset_list3, s_coeffs, 
             label=r'$\delta V4=%.2f uV$' % (round(offs_inputs[4], 2)))
    
plt.plot(offset_list3, true_seebeck, 'r--', label="True Seebeck Coefficient")

plt.fill_between(offset_list3, minus_percent, plus_percent, color='gray', 
                 alpha=.3, label="$\u00B1$%d"%(percent_error)+
                 "% Typical Experimental\nUncertainty")
plt.title("Offsets in Cold Thermocouple Only, Tref=%d K"%(T_ref_K), pad=20)
plt.xlabel('$\delta$V3 ($\mu$V)', fontsize=font_size)
plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=font_size)
plt.legend(bbox_to_anchor=(1.05,1))
invoke_plot_params("S_vs_dV3_dV4")

# hold offsets 2 and 4 constant while varying 1 and 3------------ FIGURE 3
offset_list1 = offset_list3 # offset_list1 becomes a linspace for plot
offset_list3 = main_offset_list # reset offset_list3
offs_inputs = [0, 0, 0, 0, 0]
for ind in range(len(offset_list3)):
    s_coeffs = []
    offs_inputs[3] = offset_list3[ind]
    for offset1 in offset_list1:
        offs_inputs[1] = offset1
        s_coeffs.append(seebeck_measurement(Thots_C, Tcolds_C, offs_inputs))

    plt.plot(offset_list1, s_coeffs, 
             label=r'$\delta V3=%.2f uV$' % (round(offs_inputs[3], 2)))
    
plt.plot(offset_list1, true_seebeck, 'r--', label="True Seebeck Coefficient")

plt.fill_between(offset_list1, minus_percent, plus_percent, color='gray', 
                 alpha=.3, label="$\u00B1$%d"%(percent_error)+
                 "% Typical Experimental\nUncertainty")
plt.title("Offsets in Constantan Wires Only, Tref=%d K"%(T_ref_K), pad=20)
plt.xlabel('$\delta$V1 ($\mu$V)', fontsize=font_size)
plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=font_size)
plt.legend(bbox_to_anchor=(1.05,1))
invoke_plot_params("S_vs_dV1_dV3")

# hold offsets 1 and 3 constant while varying 2 and 4------------ FIGURE 4
offset_list2 = offset_list1 # offset_list2 becomes a linspace for plot
offset_list1 = main_offset_list # reset offset_list1
offs_inputs = [0, 0, 0, 0, 0]
for ind in range(len(offset_list4)):
    s_coeffs = []
    offs_inputs[4] = offset_list4[ind]
    for offset2 in offset_list2:
        offs_inputs[2] = offset2
        s_coeffs.append(seebeck_measurement(Thots_C, Tcolds_C, offs_inputs))
    
    plt.plot(offset_list2, s_coeffs, 
             label=r'$\delta V4=%.2f uV$' % (round(offs_inputs[4], 2)))
    
plt.plot(offset_list2, true_seebeck, 'r--', label="True Seebeck Coefficient")

plt.fill_between(offset_list2, minus_percent, plus_percent, color='gray', 
                 alpha=.3, label="$\u00B1$%d"%(percent_error)+
                 "% Typical Experimental\nUncertainty")
plt.title("Offsets in Copper Wires Only, Tref=%d K"%(T_ref_K), pad=20)
plt.xlabel('$\delta$V2 ($\mu$V)', fontsize=font_size)
plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=font_size)
plt.legend(bbox_to_anchor=(1.05,1))
invoke_plot_params("S_vs_dV2_dV4")

# plot true/measured seebeck voltages vs. true/measured temperature differences
if False: # only plot if needed
    offs_inputs = [0, 50, 50, 0, 0]
    # use conversion polynomials to get delta V values
    delta_V12_true = [temp_to_voltage(temp, Tref_C) for temp in Thots_C] # mV
    delta_V34_true = [temp_to_voltage(temp, Tref_C) for temp in Tcolds_C] # mV
    # add simulated voltage offsets, convert to mV
    delta_V12_meas = [volt + offs_inputs[1]/1000 + offs_inputs[2]/1000 
                      for volt in delta_V12_true] # mV
    delta_V34_meas = [volt + offs_inputs[3]/1000 + offs_inputs[4]/1000 
                      for volt in delta_V34_true] # mV
    
    # use polynomials to return to temperatures
    offs_Thots_C = [voltage_to_temp(volt, Tref_C) for volt in delta_V12_meas]
    offs_Tcolds_C = [voltage_to_temp(volt, Tref_C) for volt in delta_V34_meas]
    meas_dT = [offs_Thots_C[ind]-offs_Tcolds_C[ind] for
          ind in range(len(offs_Thots_C))]
    
    S_Cu = round(Seebeck_Cu(T_ref_K), 3) # units: uV/K
    S_Con = round(Seebeck_constantan(T_ref_K), 3) # units: uV/K
    true_deltaV13 = [-1*(Seebeck_SRM3451(T_ref_K) - S_Con)*delta_T for
                     delta_T in dT_true] # uV
    true_deltaV24 = [-1*(Seebeck_SRM3451(T_ref_K) - S_Cu)*delta_T for
                     delta_T in dT_true] # uV
    # note: true_deltaV is in uV
    # introduce voltage offset for true_deltaV lists
    meas_deltaV13 = [volt + offs_inputs[1] + offs_inputs[3] 
                     for volt in true_deltaV13] # uV
    meas_deltaV24 = [volt + offs_inputs[2] + offs_inputs[4] 
                     for volt in true_deltaV24] # uV
    meas_deltaV = [meas_deltaV24, meas_deltaV13] # uV
    use_top_13_wires = False # choose between deltaV13 or deltaV24 for Seebeck voltage
    
    true_trend = calculate_trendline(dT_true, true_deltaV24)
    meas_trend = calculate_trendline(meas_dT, meas_deltaV[use_top_13_wires])
    
    plt.plot(dT_true, true_deltaV24, label="True Values")
    plt.plot(meas_dT, meas_deltaV[use_top_13_wires], label="Measured Values")
    plt.title('Voltage vs. Temperature, Measured and True', pad=20)
    plt.xlabel(r'True/Measured $\Delta$ Temperature (K)', fontsize=font_size)
    plt.ylabel('Seebeck Voltgae (uV)', fontsize=font_size)
    plt.legend(bbox_to_anchor=(1.05,1))
    invoke_plot_params("true_meas_dV_vs_dT")

# plot meas_deltaV13 vs. meas_deltaV24 for various offsets
if False: # only plot if needed
    offs = [0,0,0,0,0] # initialize new offsets specific to this graph
    offs_combos = [(-100,-100), (-100,0), (0,0), (100,0), (100,100)]
    S_Cu = round(Seebeck_Cu(T_ref_K), 3) # units: uV/K
    S_Con = round(Seebeck_constantan(T_ref_K), 3) # units: uV/K
    true_deltaV13 = [-1*(Seebeck_SRM3451(T_ref_K) - S_Con)*delta_T for
                     delta_T in dT_true]
    true_deltaV24 = [-1*(Seebeck_SRM3451(T_ref_K) - S_Cu)*delta_T for
                     delta_T in dT_true] #dT_true is a global variable
    # note: true_deltaV is in uV
    for tup in offs_combos:
        offs[1] = tup[0]
        offs[3] = tup[1]
        # introduce voltage offset for true_deltaV lists
        meas_deltaV13 = [volt + offs[1] + offs[3] for volt in true_deltaV13]
        meas_deltaV24 = [volt + offs[2] + offs[4] for volt in true_deltaV24]
        # both meas_deltaV are in uV
        dV_trend = calculate_trendline(meas_deltaV24, meas_deltaV13)
        plt.plot(meas_deltaV24, meas_deltaV13, 
                 label=r'$\delta V_{1}=%.2f uV$   $\delta V_{3}=%.2f uV$' % 
                 (offs[1],offs[3]) + "\nslope = %.5f"%(dV_trend['slope'])
                 + "\nintercept = %.2f"%(dV_trend['intercept']))
    plt.title("Constantan Terminal Voltage vs. Copper Terminal Voltage"+
              ", Tref=%d K"%(T_ref_K), pad=20)
    plt.xlabel('$\Delta V_{24}$ Measured ($\mu$V)', fontsize=font_size)
    plt.ylabel('$\Delta V_{13}$ Measured ($\mu$V)', fontsize=font_size)
    plt.legend(bbox_to_anchor=(1.05,1))
    invoke_plot_params("V13_vs_V24_meas")

if False: #print two thermocouple voltage values, calculate a realistic case
    offs = [0, 50, 0, 50, 0]
    seebeck_measurement(Thots_C, Tcolds_C, offs, plot=True, print_vals=True)

if False: # plot temperature offset in Delta T values for increasing Delta T
    offs_inputs = [0, 50, 50, 0, 0]
    offs_Delta_T = plot_DT_offset(offs_inputs, T_ref_K)
    plt.plot(dT_true, offs_Delta_T, 'k')
    plt.title('Offset in $\Delta T$ vs. True $\Delta T$', pad=20)
    plt.xlabel(r'True $\Delta T$ (K)', fontsize=font_size)
    plt.ylabel(r'$\delta (\Delta T)$ (K)', fontsize=font_size)
    invoke_plot_params("offs_DT_vs_DTtrue")
    
if False: # plot polynomials and show greater offset effect at higher temps
    # plot_polynomials([-200, 100], [-4, 20], Tref_C) # wide ranges
    # plot_polynomials([-200, -180], [-0.1, 0.3], Tref_C) # realistic ranges
    plot_polynomials([-200, -180], [0, 0.15], Tref_C) # realistic ranges
    


















