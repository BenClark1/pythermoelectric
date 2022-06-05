# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 21:15:09 2021

Simulated measurement of Seebeck Coefficient including affects of spurious
voltage offsets. Capable of simulating type T thermocouples with NIST SRM3451
(Bi2Te3+x) or type R thermocouples with NIST SRM3452 (Si80Ge20)

@author: benja
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import math
from decimal import Decimal
from statistics import mean
import numpy as np

import thermocouple_coefficients as tcc


# functions ----------------------------------------
def temp_to_voltage(temp_in, T_ref):
    # convert temperature (C) to EMF (mV) type T thermocouple
    if THERMOCOUPLE_TYPE == 'type T':
        # type T thermocouple: only for -270C <= temp_in <= 400C
        emf_pre = 0 # preliminary voltage before subtracting out emf_ref
        if temp_in > 0 and temp_in <= 400:
            for ind in range(len(tcc.b_pos_typeT)):
                emf_pre += tcc.b_pos_typeT[ind]*(temp_in**ind)
        elif temp_in >= -270 and temp_in <= 0:  # for negative temp inputs
            for ind in range(len(tcc.b_neg_typeT)):
                emf_pre += tcc.b_neg_typeT[ind]*(temp_in**ind)
        else:
            raise ValueError('temperature cannot be converted, out of range')

        emf_ref = 0
        if T_ref > 0 and T_ref <= 400:
            for ind in range(len(tcc.b_pos_typeT)):
                emf_ref += tcc.b_pos_typeT[ind]*(T_ref**ind)
        elif T_ref >= -270 and T_ref <= 0:  # for negative temp inputs
            for ind in range(len(tcc.b_neg_typeT)):
                emf_ref += tcc.b_neg_typeT[ind]*(T_ref**ind)
        else:
            raise ValueError('temperature cannot be converted, out of range')


    elif THERMOCOUPLE_TYPE == 'type R':
        # type R thermocouple: only for -50C <= temp_in <= 1768.1C
        emf_pre = 0 # preliminary voltage before subtracting out emf_ref
        if temp_in >= -50 and temp_in <= 1064.18:  # range 1 (r1)
            for ind in range(len(tcc.b_r1_typeR)):
                emf_pre += tcc.b_r1_typeR[ind]*(temp_in**ind)
        elif temp_in > 1064.18 and temp_in <= 1664.5:  # range 2 (r2)
            for ind in range(len(tcc.b_r2_typeR)):
                emf_pre += tcc.b_r2_typeR[ind]*(temp_in**ind)
        elif temp_in > 1664.5 and temp_in <= 1768.1:  # range 3 (r3)
            for ind in range(len(tcc.b_r3_typeR)):
               emf_pre += tcc.b_r3_typeR[ind]*(temp_in**ind)
        else:
            raise ValueError('temperature cannot be converted, out of range')

        emf_ref = 0
        if T_ref >= -50 and T_ref <= 1064.18:   # range 1 (r1)
            for ind in range(len(tcc.b_r1_typeR)):
                emf_ref += tcc.b_r1_typeR[ind]*(T_ref**ind)
        elif T_ref > 1064.18 and T_ref <= 1664.5:
            for ind in range(len(tcc.b_r2_typeR)):  # range 2 (r2)
                emf_ref += tcc.b_r2_typeR[ind]*(T_ref**ind)
        elif T_ref > 1664.5 and T_ref <= 1768.1:
            for ind in range(len(tcc.b_r3_typeR)):  # range 3 (r3)
                emf_ref += tcc.b_r3_typeR[ind]*(T_ref**ind)
        else:
            raise ValueError('temperature cannot be converted, out of range')


    else:
        raise ValueError('global constant THERMOCOUPLE_TYPE is incorrect')


    return emf_pre - emf_ref # returns emf in millivolts (mV)

def voltage_to_temp(emf_measured, T_ref):
    # convert EMF (mV) to temperature (C) for type T thermocouple
    if THERMOCOUPLE_TYPE == 'type T':
    # type T thermocouple: only for -5.603mV <= emf_measured <= 20.872mV
        temp_out = 0
        emf_ref = 0 # compute EMF for reference temperature using polynomial
        if T_ref >= 0:
            for ind in range(len(tcc.b_pos_typeT)):
                emf_ref += tcc.b_pos_typeT[ind]*(T_ref**ind)
        else:
            for ind in range(len(tcc.b_neg_typeT)):
                emf_ref += tcc.b_neg_typeT[ind]*(T_ref**ind)

        final_emf = emf_measured + emf_ref
        if final_emf >= 0:
            for ind in range(len(tcc.c_pos_typeT)):
                temp_out += tcc.c_pos_typeT[ind]*(final_emf**ind)
        else:
            for ind in range(len(tcc.c_neg_typeT)):
                temp_out += tcc.c_neg_typeT[ind]*(final_emf**ind)


    elif THERMOCOUPLE_TYPE == 'type R':
        # type R thermocouple: only for -0.226mV to 21.103mV
        temp_out = 0
        emf_ref = 0
        if T_ref >= -50 and T_ref <= 1064.18:   # range 1 (r1)
            for ind in range(len(tcc.b_r1_typeR)):
                emf_ref += tcc.b_r1_typeR[ind]*(T_ref**ind)
        elif T_ref > 1064.18 and T_ref <= 1664.5:
            for ind in range(len(tcc.b_r2_typeR)):  # range 2 (r2)
                emf_ref += tcc.b_r2_typeR[ind]*(T_ref**ind)
        elif T_ref > 1664.5 and T_ref <= 1768.1:
            for ind in range(len(tcc.b_r3_typeR)):  # range 3 (r3)
                emf_ref += tcc.b_r3_typeR[ind]*(T_ref**ind)
        else:
            raise ValueError(
                'reference temperature out of range for type R thermocouple')

        final_emf = emf_measured + emf_ref

        if final_emf >= -0.226 and final_emf <= 1.923:  # range 1 (r1)
            for ind in range(len(tcc.c_r1_typeR)):
                temp_out += tcc.c_r1_typeR[ind]*(final_emf**ind)
        elif final_emf > 1.923 and final_emf <= 13.228:  # range 2 (r2)
            for ind in range(len(tcc.c_r2_typeR)):
                temp_out += tcc.c_r2_typeR[ind]*(final_emf**ind)
        # 13.228, 11.361: overlapping ranges as per NIST polynomial definition
        elif final_emf > 11.361 and final_emf <= 19.739:  # range 3 (r3)
            for ind in range(len(tcc.c_r3_typeR)):
                temp_out += tcc.c_r3_typeR[ind]*(final_emf**ind)
        elif final_emf > 19.739 and final_emf <= 21.103:  # range 4 (r4)
            for ind in range(len(tcc.c_r4_typeR)):
                temp_out += tcc.c_r4_typeR[ind]*(final_emf**ind)
        else:
            if final_emf < -0.226:
                raise ValueError(
                    'voltage cannot be converted, less than -0.226mV')
            else:
                raise ValueError(
                    'voltage cannot be converted, more than 21.103mV')
    else:
        raise ValueError('global constant THERMOCOUPLE_TYPE incorrect')

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


def seebeck_Cu(T):
# T must be in kelvin: Seebeck coeff [uV/K]
    if 73.15 < T and T < 673.15:
        exp_term = math.exp(-1*T/93)
        rational_fraction = 0.442/(1+(T/172.4)**3)
        return 0.041*T*( exp_term + 0.123 -  rational_fraction) + 0.804
    else:
        raise ValueError("can't get Copper Seebeck, temperature out of range")

def seebeck_SRM3451(T):
    # Standard Reference Material (SRM) 3451 = Bi2Te3+x
    # NIST provided Seebeck coefficient for SRM3451
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

def seebeck_SRM3452_SiGe(T):  # Si80Ge20 Seebeck coefficient
    # T must be in KELVIN: Seebeck coeff [uV/K]
    A = 295
    S_A = 1.16246764e2  # uV / K
    a = 2.343158e-1
    b = -8.781594e-5
    term2 = a * (T - A)
    term3 = b * (T - A)**2
    return S_A + term2 + term3

def seebeck_constantan(T): # units of T: K
    # get the Seebeck coefficent of Constantan
    # first get Seebeck coefficient for Copper/Constantan T-type thermocouple
    if 73.15 < T and T < 673.15: # MUST BE FOR 73.15 K < T < 673.15 K
        S_Cu_Con = 4.37184 + 0.1676*T - (1.84371e-4)*T**2 \
            + (1.2244e-7)*T**3 - (4.47618e-11)*T**4
        # subtract Cu/Con Seebeck value from copper Seebeck to get Con Seebeck
        # S_Cu_Con is the type T thermocouple Seebeck coefficient
        return seebeck_Cu(T) - S_Cu_Con
    else:
        raise ValueError("can't get Constantan Seebeck, temperature out of range")


def seebeck_platinum(T): # units of T: K
    # get the Seebeck coefficient of Platinum
    if 70 < T and T < 1500:
        bracket_term = (math.exp(-1*T/88) - 0.0786 + 0.43 / (1 + (T/84.3)**4))
        return 0.186 * T * bracket_term - 2.57
    else:
        raise ValueError("can't get platinum Seebeck, temperature out of range")

def seebeck_measurement(
        Thots_C, Tcolds_C, offs, tref_K,
        plot=False, print_vals=False) -> tuple[float, bool]:
    # the tuple[float, bool] is just documentation of the return type
    # indicate if polynomials are acting near the transition values
    tref_C = kelvin_to_celsius(tref_K)
    # offs: a list of 5 offsets, first term must be zero
        # index of offs corresponds to offset location e+
    # use conversion polynomials to get delta V values
    delta_V12_true = [temp_to_voltage(temp, tref_C) for temp in Thots_C] # mV
    delta_V34_true = [temp_to_voltage(temp, tref_C) for temp in Tcolds_C] # mV

    in_transition = False
    if enable_check_transitions:  # transitions only supported for type R
        transition_temps = [-50, 1064.18, 1664.5, 1768.1]
        # transition_temps = [-50, 127, 1664.5, 1768.1]  # temporary test
        for trans_temp in transition_temps:
            if Tcolds_C[0] < trans_temp and Tcolds_C[-1] > trans_temp:
                print("found temperature transition")
                in_transition = True
            if Thots_C[0] < trans_temp and Thots_C[-1] > trans_temp:
                in_transition = True
                print("found temperature transition")

    # add simulated voltage offsets, convert to mV
    delta_V12_meas = [volt + offs[1]/1000 + offs[2]/1000 for volt in delta_V12_true]
    delta_V34_meas = [volt + offs[3]/1000 + offs[4]/1000 for volt in delta_V34_true]
    if print_vals: # print two thermocouple voltages to see realistic example
        print("Printing voltages (mv)   tref=%.1fK   thermocouple: %s...\n" %
              (TREF_K, THERMOCOUPLE_TYPE))
        print("true hot voltage low power: ", round(delta_V12_true[1], 5))
        # print("offset hot voltage low power: ", round(delta_V12_meas[1], 5))
        print("true cold voltage low power: ", round(delta_V34_true[1], 5))
        # print("offset cold voltage low power: ", round(delta_V34_meas[1], 5))
        print("true hot voltage high power: ", round(delta_V12_true[9], 5))
        # print("offset hot voltage high power: ", round(delta_V12_meas[8], 5))
        print("true cold voltage high power: ", round(delta_V34_true[9], 5))
        # print("offset cold voltage high power: ", round(delta_V34_meas[8], 5))
    # use polynomials to return to temperatures
    offs_Thots_C = [voltage_to_temp(volt, tref_C) for volt in delta_V12_meas]
    offs_Tcolds_C = [voltage_to_temp(volt, tref_C) for volt in delta_V34_meas]

    if enable_check_transitions:
        transition_volts = [-0.226, 1.923, 13.228, 19.739, 21.103]
        for trans_volt in transition_volts:
            if delta_V12_meas[0] < trans_volt and delta_V12_meas[-1] > trans_volt:
                print("found voltage transition")
                in_transition = True
            if delta_V34_meas[0] < trans_volt and delta_V34_meas[-1] > trans_volt:
                in_transition = True
                print("found voltage transition")

    # meas_dT is the same in both Kelvin and Celsius
    meas_dT = [offs_Thots_C[ind]-offs_Tcolds_C[ind] for
              ind in range(len(offs_Thots_C))]

    S_Cu = seebeck_Cu(tref_K) # units: uV/K
    S_Con = seebeck_constantan(tref_K) # units: uV/K
    S_Pt = seebeck_platinum(tref_K) # units: uV/K
    if THERMOCOUPLE_TYPE == 'type T':  # true_deltaV13 is different for type T and R
        true_deltaV13 = [-1*(seebeck_SRM3451(tref_K) - S_Con)*delta_T for
                         delta_T in dT_true]
        true_deltaV24 = [-1*(seebeck_SRM3451(tref_K) - S_Cu)*delta_T for
                         delta_T in dT_true] #dT_true is a global variable
        # note: true_deltaV is in uV
        # introduce voltage offset for true_deltaV lists
        meas_deltaV13 = [volt + offs[1] + offs[3] for volt in true_deltaV13] # uV
        meas_deltaV24 = [volt + offs[2] + offs[4] for volt in true_deltaV24] # uV
        # choose between deltaV13 or deltaV24 for Seebeck voltage
        meas_deltaV = [meas_deltaV24, meas_deltaV13] # uV

        use_top_13_wires = False

    elif THERMOCOUPLE_TYPE == 'type R':
        true_deltaV13 = [-1*(seebeck_SRM3452_SiGe(tref_K) - S_Pt)*delta_T for
                         delta_T in dT_true]
        # note: true_deltaV is in uV
        # introduce voltage offset for true_deltaV lists
        meas_deltaV13 = [volt + offs[1] + offs[3] for volt in true_deltaV13] # uV
        # meas_deltaV24 = [volt + offs[2] + offs[4] for volt in true_deltaV24] # uV
        meas_deltaV = [None, meas_deltaV13] # uV
        # choose between deltaV13 or deltaV24 for Seebeck voltage
        # NOTE: always use 13 for type R because we don't have S for Pt-13%Rh
        use_top_13_wires = True
    else:
        raise ValueError('global constant THERMOCOUPLE_TYPE is incorrect')

    # get a dictionary with slope, intercept, and trendline y values
    trend_info = calculate_trendline(meas_dT, meas_deltaV[use_top_13_wires])

    if plot:
        plt.plot(meas_dT, meas_deltaV[use_top_13_wires], 'b.')
        plt.plot(meas_dT, trend_info['trendline'], 'b')
        # plt.title('Thermoelectric Votlage Produced by Seebeck Effect in Bi₂Te₃₊ₓ',
        #           pad=20)
        plt.title('Thermoelectric Votlage Produced by Seebeck Effect',
                  pad=20)
        plt.xlabel('Measured $\Delta$ Temperature (K)', fontsize=font_size)
        plt.ylabel('Measured Seebeck Voltage (uV)', fontsize=font_size)
        invoke_plot_params("Seebeck_Effect_for_Slope")

    if THERMOCOUPLE_TYPE == 'type T':
        S_sample = -1*trend_info['slope'] + [S_Cu, S_Con][use_top_13_wires]
    elif THERMOCOUPLE_TYPE == 'type R':
        S_sample = -1*trend_info['slope'] + [None, S_Pt][use_top_13_wires]
    # print("\nFinal Seebeck Coefficient of the Sample: ")
    # print(round(S_sample, 9))

    # print("\nDifferences between original and new temps (C): ")
    # print("hot: ")
    # print([Thots_C[i]-offs_Thots_C[i] for i in range(len(offs_Thots_C))])
    # print("cold: ")
    # print([Tcolds_C[i]-offs_Tcolds_C[i] for i in range(len(offs_Tcolds_C))])

    return (S_sample, in_transition)


def invoke_plot_params(filename="unnamed_1200dpi"):
    plt.tick_params(axis='both',which='major', direction = 'in',
                    top = 1, right = 1, length=6,width=1,labelsize=font_size)
    plt.tick_params(axis='both',which='minor', direction = 'in',
                    top = 1, right = 1, length=2,width=1,labelsize=font_size)
    plt.grid()
    if SAVEFIG:
        # plt.savefig("plot_img_1200dpi/"+filename, dpi=1200, bbox_inches='tight')
        plt.savefig("post_APS_plots/"+filename, dpi=1200, bbox_inches='tight')
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


def plot_polynomials(temp_range, volt_range, plot_T_to_V=False):
    # temp_range units: C   volt_range units: mV
    if plot_T_to_V:
        temps_in = np.linspace(temp_range[0], temp_range[1], NUM_POINTS) # units: C
        volts_out = [temp_to_voltage(temp, TREF_C) for temp in temps_in]

        plt.plot(temps_in, volts_out, 'k')
        plt.title("Temperature to Voltage Conversion Polynomial, " +
                  "Tref=%d K, Thermocouple %s" %
                  (TREF_K, THERMOCOUPLE_TYPE), pad=20)
        # plt.subtitle('test')
        plt.xlabel('Temperature (C)', fontsize=font_size)
        plt.ylabel('Voltage (mV)', fontsize=font_size)
        invoke_plot_params("temp_to_voltage")

    volts_in = np.linspace(volt_range[0], volt_range[1], NUM_POINTS) #units: mV
    temps_out = [voltage_to_temp(volt, TREF_C) for volt in volts_in]

    # LP stands for low power
    # volt_diff = 0.05 # mV  # this represents a spurious voltage offset
    volt_diff = 0.02 # mV  # this represents a spurious voltage offset
    cold_start_LP = 0.00348 # mV  # FOR TREF=80K, TYPE T, 1mW
    # cold_start_LP = 0.00697 # mV  # FOR TREF=80K, TYPE T, 2 mW
    # cold_start_LP = 0.01396 # mV  # FOR TREF=80K, TYPE T, 4 mW
    cold_volts_LP = [cold_start_LP, cold_start_LP + -1*volt_diff]
    hot_start_LP = 0.00697 # mV  # FOR TREF=80K, TYPE T, 1mW
    # hot_start_LP = 0.01396 # mV  # FOR TREF=80K, TYPE T, 2 mW
    # hot_start_LP = 0.02801 # mV  # FOR TREF=80K, TYPE T, 4 mW
    hot_volts_LP = [hot_start_LP, hot_start_LP + volt_diff]
    cold_true_LP = voltage_to_temp(cold_volts_LP[0], TREF_C)
    cold_meas_LP = voltage_to_temp(cold_volts_LP[1], TREF_C)
    hot_true_LP = voltage_to_temp(hot_volts_LP[0], TREF_C)
    hot_meas_LP = voltage_to_temp(hot_volts_LP[1], TREF_C)
    cold_temps_LP = [cold_true_LP, cold_meas_LP]
    hot_temps_LP = [hot_true_LP, hot_meas_LP] # 1 is "true" 2 is "measured" after offset

    # cold_start = 0.02801 # mV  # FOR TREF=80K, TYPE T, 8 mW
    cold_start = 0.03153 # mV  # FOR TREF=80K, TYPE T, 9 mW
    cold_volts = [cold_start, cold_start + -1*volt_diff]
    # hot_start = 0.05639 # mV  # FOR TREF=80K, TYPE T, 8 mW
    hot_start = 0.06355 # mV  # FOR TREF=80K, TYPE T, 9 mW
    hot_volts = [hot_start, hot_start + volt_diff]
    cold_true = voltage_to_temp(cold_volts[0], TREF_C)
    cold_meas = voltage_to_temp(cold_volts[1], TREF_C)
    hot_true = voltage_to_temp(hot_volts[0], TREF_C)
    hot_meas = voltage_to_temp(hot_volts[1], TREF_C)
    cold_temps = [cold_true, cold_meas]
    hot_temps = [hot_true, hot_meas]  # 1 is "true" 2 is "measured" after offset

    plt.plot(volts_in, temps_out, 'k')  # ploting polynomial curve

    # plt.plot(cold_volts_LP, cold_temps_LP, 'o', color='green',
    #           label='Temp diff: %.3f C' % (cold_meas_LP - cold_true_LP))
    # plt.plot(hot_volts_LP, hot_temps_LP, 'o', color='orange',
    #           label='Temp diff: %.3f C' % (hot_meas_LP-hot_true_LP))
    plt.plot(cold_volts_LP, cold_temps_LP, 'bo', markerfacecolor='none',
              label='Temp diff: %.3f C' % (cold_meas_LP - cold_true_LP))
    plt.plot(hot_volts_LP, hot_temps_LP, 'ro', markerfacecolor = 'none',
              label='Temp diff: %.3f C' % (hot_meas_LP-hot_true_LP))

    true_diff_LP = hot_temps_LP[0] - cold_temps_LP[0]
    meas_diff_LP = hot_temps_LP[1] - cold_temps_LP[1]
    # err_in_DT is the error in delta T
    err_in_DT_LP = meas_diff_LP - true_diff_LP
    print("Low power measurements (C): ")
    print("\tTrue temperature difference: %.4f" % true_diff_LP)
    print("\tMeasured temperature difference: %.4f" % meas_diff_LP)
    print("\tError in DT: %.4f" % err_in_DT_LP)

    plt.title("Voltage to Temperature Conversion Polynomial, " +
              "Tref=%d K, Thermocouple %s" %
              (TREF_K, THERMOCOUPLE_TYPE), pad=20)
    plt.xlabel('Voltage (mV)', fontsize=font_size)
    plt.ylabel('Temperature (C)', fontsize=font_size)
    invoke_plot_params("voltage_to_temp_LP")

    plt.plot(volts_in, temps_out, 'k')  # ploting polynomial curve

    plt.plot(cold_volts, cold_temps, 'bo',
              label='Temp diff: %.3f C' % (cold_meas-cold_true))
    plt.plot(hot_volts, hot_temps, 'ro',
              label='Temp diff: %.3f C' % (hot_meas-hot_true))

    true_diff = hot_temps[0] - cold_temps[0]
    meas_diff = hot_temps[1] - cold_temps[1]
    err_in_DT = meas_diff - true_diff
    print("High power measurements (C): ")
    print("\tTrue temperature difference: %.4f" % true_diff)
    print("\tMeasured temperature difference: %.4f" % meas_diff)
    print("\tError in DT: %.4f" % err_in_DT)

    # plt.axvline(0, -195, -180, color='r') # for plotting horizontal lines
    # test = plt.axvline(x=.2, ymin=-195, ymax=-180)
    # plt.plot(test)
    plt.title("Voltage to Temperature Conversion Polynomial, " +
              "Tref=%d K, Thermocouple %s" %
              (TREF_K, THERMOCOUPLE_TYPE), pad=20)
    plt.xlabel('Voltage (mV)', fontsize=font_size)
    plt.ylabel('Temperature (C)', fontsize=font_size)
    # plt.legend(bbox_to_anchor=(1.05,1))
    invoke_plot_params("voltage_to_temp")

# main script ---------------------------------------

# CONTROL PANEL: configure these settings before running the code
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# reference temperature at the base of the sample
# TREF_K = 80 # units: K
# TREF_K = 293 # units: K
# TREF_K = 270 # units: K
# TREF_K = 304 # units: K
# TREF_K = 400 # units: K
TREF_K = 500 # units: K
# TREF_K = 520 # units: K
# TREF_K = 670 # units: K  # max value until temp out of range for Cu seebeck
TREF_C = kelvin_to_celsius(TREF_K) # reference temperature in celsius

# create offsets in uV
main_offset_list = [-200, -100, -50, 0, 50, 100, 200]
# main_offset_list = [-400, -100, -50, 0, 50, 100, 400]
offset_list1 = main_offset_list
offset_list2 = main_offset_list
offset_list3 = main_offset_list
offset_list4 = main_offset_list
ind_zero = len(offset_list1)//2  # get index of zero offset (center of list)

font_size = 16
# all caps - constants should not be changed in runtime
PERCENT_ERROR = 5  # %  5% is typical measurement uncertainty for seebeck coef
NUM_POINTS = 301  # choose the resolution of horizontal axis for plots
# NUM_POINTS = 21
THERMOCOUPLE_TYPES = ('type T', 'type R')  # do not change this line
# never change THERMOCOUPLE_TYPE outside of this line
THERMOCOUPLE_TYPE = THERMOCOUPLE_TYPES[1]  # change the index here when needed
# SAVEFIG will save any figures as png files when invoke_plot_params() is called
SAVEFIG = False

enable_seebeck_vs_offsets_plots = True  # enable plots as needed
enable_seebeck_volts_vs_temp_diff = False
enable_measdV13_vs_measdV24 = False
enable_print_TC_volts = False
enable_dDT_vs_trueDT = False
enable_plot_polynomials = False
enable_percent_error_plot = True
enable_plot_blip = False
enable_check_transitions = False # transitions only supported for type R
enable_plot_dashed_DV_vs_DT = False
enable_examine_blip_region = False

print("Running simulation code...")
print("Reference temperature: %.2f K = %.2f C" % (TREF_K, TREF_C))
print("Thermocouple type: ", THERMOCOUPLE_TYPE, '\n')
# ------------------------------------------------------------------------
# end control panel ------------------------------------------------------

# dQ/dT_true = P = kA(dT_true/dx) = power delivered to the sample
powers = [0, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, \
          6e-3, 7e-3, 8e-3, 9e-3, 10e-3] # units: W eventually test 100 values here
# powers = np.linspace(0, 10e-3, 101) # units: W
kappa = 2.0 # units: W/(m*K)  thermal conductivity kappa for Bi2Te3
area = 0.002**2 # units: m^2  2mm x 2mm cross sectional area
# location of hot and cold thermocoulpe probe points
x_cold = 0.001667 # 1.67mm, or sample length / 3
dx = 0.001667 # 1.67mm  length difference between each thermocouple probe point
x_hot = x_cold + dx # 2 * 1.67mm
# difference in temperature between each thermocouple (delta T)
dT_true = [(pwr * dx)/(kappa * area) for pwr in powers] #units: K
# derivative of T with respect to x is
    # dT_true/dx (slope)  T(x) = (dT_true/dx)x + TREF_K

# get hot and cold temperatures and convert to celsius
Thots = [(delta_T/dx)*x_hot + TREF_K for delta_T in dT_true]
Tcolds = [(delta_T/dx)*x_cold + TREF_K for delta_T in dT_true]
Thots_C = [kelvin_to_celsius(temp) for temp in Thots]
Tcolds_C = [kelvin_to_celsius(temp) for temp in Tcolds]

if THERMOCOUPLE_TYPE == 'type T':  # keep track of material of each wire number
    material_string_13 = "Constantan"  # 1 and 3 are the negative leads
    material_string_24 = "Copper"  # 2 and 4 are the positive leads
    true_seebeck = [seebeck_SRM3451(TREF_K)] * NUM_POINTS
elif THERMOCOUPLE_TYPE == 'type R':
    material_string_13 = "Platinum"  # these are negative
    material_string_24 = "Pt-13% Rh"  # these are positive
    true_seebeck = [seebeck_SRM3452_SiGe(TREF_K)] * NUM_POINTS
else:
    raise ValueError('global constant THERMOCOUPLE_TYPE is incorrect')

minus_percent = [S*(1-(PERCENT_ERROR/100)) for S in true_seebeck]
plus_percent = [S*(1+(PERCENT_ERROR/100)) for S in true_seebeck]

offs_inputs = [0, 0, 0, 0, 0]

round_and_print("Applying the following powers (W) to the heater: ", powers, 7)

# ------------------------------------------------------------------------
if enable_seebeck_vs_offsets_plots:
    # hold offsets 3 and 4 constant while varying 1 and 2 ------------ FIGURE 1
    offset_list1 = np.linspace(main_offset_list[0],main_offset_list[-1],NUM_POINTS)
    for ind in range(len(offset_list2)):
        s_coeffs = []
        offs_inputs[2] = offset_list2[ind]
        for offset1 in offset_list1:
            offs_inputs[1] = offset1
            s_coeff, in_trans = seebeck_measurement(
                Thots_C, Tcolds_C, offs_inputs, TREF_K)
            s_coeffs.append(s_coeff)
            if in_trans:
                plt.plot(offset1, s_coeff, 'kx')

        plt.plot(offset_list1, s_coeffs,
                 label=r'$\delta V2=%.2f uV$' % (round(offs_inputs[2], 2)))


    plt.plot(offset_list1, true_seebeck, 'r--', label="True Seebeck Coefficient")
    # plt.plot(offset_list1, minus_percent, 'k:',
    #          label='%d'%(PERCENT_ERROR)+"% Error Bounds")
    # plt.plot(offset_list1, plus_percent, 'k:')
    plt.fill_between(offset_list1, minus_percent, plus_percent, color='gray',
                     alpha=.3, label="$\u00B1$%d"%(PERCENT_ERROR)+
                     "% Typical Experimental\nUncertainty")
    plt.title("Offsets in Hot Thermocouple Only, Tref=%d K"%(TREF_K), pad=20)
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
            # calculate the Seebeck coefficient and determine if it was
            # calculated over a transition boundary:
            s_coeff, in_trans = seebeck_measurement(
                Thots_C, Tcolds_C, offs_inputs, TREF_K)
            s_coeffs.append(s_coeff)
            if in_trans:
                plt.plot(offset3, s_coeff, 'kx')

        plt.plot(offset_list3, s_coeffs,
                 label=r'$\delta V4=%.2f uV$' % (round(offs_inputs[4], 2)))

    plt.plot(offset_list3, true_seebeck, 'r--', label="True Seebeck Coefficient")

    plt.fill_between(offset_list3, minus_percent, plus_percent, color='gray',
                     alpha=.3, label="$\u00B1$%d"%(PERCENT_ERROR)+
                     "% Typical Experimental\nUncertainty")
    plt.title("Offsets in Cold Thermocouple Only, Tref=%d K"%(TREF_K), pad=20)
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
            s_coeff, in_trans = seebeck_measurement(
                Thots_C, Tcolds_C, offs_inputs, TREF_K)
            s_coeffs.append(s_coeff)
            if in_trans:
                plt.plot(offset1, s_coeff, 'kx')

        plt.plot(offset_list1, s_coeffs,
                 label=r'$\delta V3=%.2f uV$' % (round(offs_inputs[3], 2)))

    plt.plot(offset_list1, true_seebeck, 'r--', label="True Seebeck Coefficient")

    plt.fill_between(offset_list1, minus_percent, plus_percent, color='gray',
                     alpha=.3, label="$\u00B1$%d"%(PERCENT_ERROR)+
                     "% Typical Experimental\nUncertainty")
    plt.title("Offsets in %s Wires Only, Tref=%d K" %
              (material_string_13, TREF_K), pad=20)
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
            s_coeff, in_trans = seebeck_measurement(
                Thots_C, Tcolds_C, offs_inputs, TREF_K)
            s_coeffs.append(s_coeff)
            if in_trans:
                plt.plot(offset2, s_coeff, 'kx')

        plt.plot(offset_list2, s_coeffs,
                 label=r'$\delta V4=%.2f uV$' % (round(offs_inputs[4], 2)))

    plt.plot(offset_list2, true_seebeck, 'r--', label="True Seebeck Coefficient")

    plt.fill_between(offset_list2, minus_percent, plus_percent, color='gray',
                     alpha=.3, label="$\u00B1$%d"%(PERCENT_ERROR)+
                     "% Typical Experimental\nUncertainty")
    plt.title("Offsets in %s Wires Only, Tref=%d K" %
              (material_string_24, TREF_K), pad=20)
    plt.xlabel('$\delta$V2 ($\mu$V)', fontsize=font_size)
    plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=font_size)
    plt.legend(bbox_to_anchor=(1.05,1))
    invoke_plot_params("S_vs_dV2_dV4")

# ------------------------------------------------------------------------
# plot true/measured seebeck voltages vs. true/measured temperature differences
if enable_seebeck_volts_vs_temp_diff: # only plot if needed
    if THERMOCOUPLE_TYPE != 'type T':
        raise ValueError("this section of code only works for type T")
    offs_inputs = [0, 50, 50, 0, 0]
    # use conversion polynomials to get delta V values
    delta_V12_true = [temp_to_voltage(temp, TREF_C) for temp in Thots_C] # mV
    delta_V34_true = [temp_to_voltage(temp, TREF_C) for temp in Tcolds_C] # mV
    # add simulated voltage offsets, convert to mV
    delta_V12_meas = [volt + offs_inputs[1]/1000 + offs_inputs[2]/1000
                      for volt in delta_V12_true] # mV
    delta_V34_meas = [volt + offs_inputs[3]/1000 + offs_inputs[4]/1000
                      for volt in delta_V34_true] # mV

    # use polynomials to return to temperatures
    offs_Thots_C = [voltage_to_temp(volt, TREF_C) for volt in delta_V12_meas]
    offs_Tcolds_C = [voltage_to_temp(volt, TREF_C) for volt in delta_V34_meas]
    meas_dT = [offs_Thots_C[ind]-offs_Tcolds_C[ind] for
          ind in range(len(offs_Thots_C))]

    S_Cu = round(seebeck_Cu(TREF_K), 3) # units: uV/K
    S_Con = round(seebeck_constantan(TREF_K), 3) # units: uV/K
    true_deltaV13 = [-1*(seebeck_SRM3451(TREF_K) - S_Con)*delta_T for
                     delta_T in dT_true] # uV
    true_deltaV24 = [-1*(seebeck_SRM3451(TREF_K) - S_Cu)*delta_T for
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

# ------------------------------------------------------------------------
# plot meas_deltaV13 vs. meas_deltaV24 for various offsets
if enable_measdV13_vs_measdV24: # only plot if needed
    if THERMOCOUPLE_TYPE != 'type T':
        raise ValueError("this section of code only works for type T")
    offs = [0,0,0,0,0] # initialize new offsets specific to this graph
    offs_combos = [(-100,-100), (-100,0), (0,0), (100,0), (100,100)]
    S_Cu = round(seebeck_Cu(TREF_K), 3) # units: uV/K
    S_Con = round(seebeck_constantan(TREF_K), 3) # units: uV/K
    true_deltaV13 = [-1*(seebeck_SRM3451(TREF_K) - S_Con)*delta_T for
                     delta_T in dT_true]
    true_deltaV24 = [-1*(seebeck_SRM3451(TREF_K) - S_Cu)*delta_T for
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
              ", Tref=%d K"%(TREF_K), pad=20)
    plt.xlabel('$\Delta V_{24}$ Measured ($\mu$V)', fontsize=font_size)
    plt.ylabel('$\Delta V_{13}$ Measured ($\mu$V)', fontsize=font_size)
    plt.legend(bbox_to_anchor=(1.05,1))
    invoke_plot_params("V13_vs_V24_meas")

# ------------------------------------------------------------------------
if enable_print_TC_volts: #print two thermocouple voltage values, calculate a realistic case
    offs = [0, 50, 0, 50, 0]
    seebeck_measurement(
        Thots_C, Tcolds_C, offs, TREF_K, print_vals=True)[0]

# ------------------------------------------------------------------------
if enable_dDT_vs_trueDT: # plot temperature offset in Delta T values for increasing Delta T
    offs_inputs = [0, 50, 50, 0, 0]
    offs_Delta_T = plot_DT_offset(offs_inputs, TREF_K)
    plt.plot(dT_true, offs_Delta_T, 'k')
    plt.title('Offset in $\Delta T$ vs. True $\Delta T$', pad=20)
    plt.xlabel(r'True $\Delta T$ (K)', fontsize=font_size)
    plt.ylabel(r'$\delta (\Delta T)$ (K)', fontsize=font_size)
    invoke_plot_params("offs_DT_vs_DTtrue")

# ------------------------------------------------------------------------
if enable_plot_polynomials: # plot polynomials and show greater offset effect at higher temps
    # plot_polynomials([-200, 100], [-4, 20]) # wide ranges
    # plot_polynomials([-200, -180], [-0.1, 0.3]) # realistic ranges
    if THERMOCOUPLE_TYPE == 'type T':
        plot_polynomials([-200, -180], [-0.02, 0.10]) # realistic ranges
    elif THERMOCOUPLE_TYPE == 'type R':
        # plot_polynomials([-40, 200], [-0.226, 21.103])
        # plot_polynomials([-40, 200], [-0.226, 0.18])
        plot_polynomials([-40, 200], [-0.02, 0.10])
        # -0.226mV to 21.103mV

# ------------------------------------------------------------------------
if enable_percent_error_plot: # plot percent error vs T for type T thermocouple
    # Tref_var is a linear space of various reference temperatures
    if THERMOCOUPLE_TYPE == 'type T':
        temp_range = [80, 400] # units: K # NOTE: reusing variable names
        Tref_var = np.linspace(temp_range[0], temp_range[1], NUM_POINTS)
        S_true = [seebeck_SRM3451(tref) for tref in Tref_var]
    elif THERMOCOUPLE_TYPE == 'type R':
        temp_range = [350, 650]
        Tref_var = np.linspace(temp_range[0], temp_range[1], NUM_POINTS)
        S_true = [seebeck_SRM3452_SiGe(tref) for tref in Tref_var]
    else:
        raise ValueError('global constant THERMOCOUPLE_TYPE is incorrect')

    # offs_var_ref contains 3 offset scenarios
    # offs_var_ref = [[0,0,0,0,0], [0,50,50,0,0], [0,-50,-50,0,0], [0,-140,-100,0,0]]
    offs_var_ref = [[0,0,0,0,0], [0,50,50,0,0], [0,-50,-50,0,0]]
    for ind in range(len(offs_var_ref)):
        S_meas = []
        for tref in Tref_var: # inefficient: nested for loops
            # thots_var is varying hot temperatures (T ref is varying)
            thots_var = [(delta_T/dx)*x_hot + tref for delta_T in dT_true]
            tcolds_var = [(delta_T/dx)*x_cold + tref for delta_T in dT_true]
            thots_C_var = [kelvin_to_celsius(temp) for temp in thots_var]
            tcolds_C_var = [kelvin_to_celsius(temp) for temp in tcolds_var]

            s_coeff, in_trans = seebeck_measurement(
                thots_C_var, tcolds_C_var, offs_var_ref[ind], tref)

            S_meas.append(s_coeff)
            if THERMOCOUPLE_TYPE == 'type T':
                if in_trans:
                    plt.plot(tref,
                             ((s_coeff - seebeck_SRM3451(tref)) * 100 /
                              seebeck_SRM3451(tref)),
                             'kx')
            elif THERMOCOUPLE_TYPE == 'type R':
                if in_trans:
                    plt.plot(tref,
                             ((s_coeff - seebeck_SRM3452_SiGe(tref)) * 100 /
                              seebeck_SRM3452_SiGe(tref)),
                             'kx')

            # S_meas.append(seebeck_measurement(thots_C_var, tcolds_C_var,
            #                                   offs_var_ref[ind], tref)[0])

        err_curve = [
            ((S_meas[ind] - S_true[ind]) * 100 / S_true[ind]) for
            ind in range(len(S_true))
            ]
        # positive percent error corresponds to "overshot" S value

        plt.plot(Tref_var, err_curve,
                 label='$\delta$V1 %d$\mu$V \n$\delta$V2 = %d$\mu$V' %
                 (offs_var_ref[ind][1], offs_var_ref[ind][2]))
        # '$\delta$V1 ($\mu$V)'

    plt.ylim([-10, 10])
    plt.title("Temperature Dependence of Seebeck Measurement Error\n" +
              "Thermocouple: %s" % THERMOCOUPLE_TYPE, pad=20)
    plt.xlabel('Reference Temperature (K)', fontsize=font_size)
    plt.ylabel('Percent Error (%)', fontsize=font_size)
    # plt.legend(bbox_to_anchor=(1.45, 1))
    plt.legend()
    invoke_plot_params("percent_error_vs_temp")

if enable_plot_blip:  # plot zoomed in Seebeck vs. offset graph to show blip
    # showing data blip for dV1 = 0, dV2 = 0, dV3 = variable, dV4 = 200uV
    blip_bounds = [-20, 20]
    offset_list3 = np.linspace(blip_bounds[0], blip_bounds[-1], NUM_POINTS)
    # for ind in range(len(offset_list2)):
    s_coeffs = []
    offs_inputs = [0, 0, 0, 0, 200]
    for offset3 in offset_list3:
        offs_inputs[3] = offset3
        s_coeffs.append(
            seebeck_measurement(Thots_C, Tcolds_C, offs_inputs, TREF_K)[0])

    plt.plot(offset_list3, s_coeffs, 'm')

    plt.title("$\delta V4=200 uV$, Tref=%d K" % (TREF_K), pad=20)
    plt.xlabel('$\delta$V3 ($\mu$V)', fontsize=font_size)
    plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=font_size)
    # plt.legend(bbox_to_anchor=(1.05,1))
    # plt.autoscale(enable=False, axis='y')

    # create a list of dV3 offsets to examine points on the blip
    blip_markers = [-18, -6, -3, 5, 10, 13.45, 18]
    for offs_mark in blip_markers:
        offs_inputs = [0, 0, 0, offs_mark, 200]
        S_result = seebeck_measurement(Thots_C, Tcolds_C, offs_inputs, TREF_K)[0]
        # plot single points for each marker
        plt.plot(offs_mark, S_result, 'r*',
                 label="$\delta$ V3=%.2f $\mu$V \nS = %.3f $\mu$V/K" %
                 (offs_mark, S_result))
    plt.legend(bbox_to_anchor=(1.05,1))
    invoke_plot_params("Single_S_zoomed_in")

    for offs_mark in blip_markers:
        offs_inputs = [0, 0, 0, offs_mark, 200]
        seebeck_measurement(Thots_C, Tcolds_C, offs_inputs, TREF_K, plot=True)[0]

if enable_plot_dashed_DV_vs_DT:
    # offs_inputs = [0, 0, 0, 0, 0]  # units: uV
    offs_inputs = [0, 20, 0, -20, 0]  # units: uV
    # replicate seebeck_measurement() but for a single plot of DV vs DT
    # use conversion polynomials to get delta V values
    delta_V12_true = [temp_to_voltage(temp, TREF_C) for temp in Thots_C] # mV
    delta_V34_true = [temp_to_voltage(temp, TREF_C) for temp in Tcolds_C] # mV

    # add simulated voltage offsets, convert to mV
    delta_V12_meas = [volt + offs_inputs[1]/1000 + offs_inputs[2]/1000 for
                      volt in delta_V12_true]
    delta_V34_meas = [volt + offs_inputs[3]/1000 + offs_inputs[4]/1000 for
                      volt in delta_V34_true]

    # use polynomials to return to temperatures
    offs_Thots_C = [voltage_to_temp(volt, TREF_C) for volt in delta_V12_meas]
    offs_Tcolds_C = [voltage_to_temp(volt, TREF_C) for volt in delta_V34_meas]

    # meas_dT is the same in both Kelvin and Celsius
    meas_dT = [offs_Thots_C[ind]-offs_Tcolds_C[ind] for
              ind in range(len(offs_Thots_C))]

    S_Cu = seebeck_Cu(TREF_K) # units: uV/K
    S_Con = seebeck_constantan(TREF_K) # units: uV/K
    S_Pt = seebeck_platinum(TREF_K) # units: uV/K
    if THERMOCOUPLE_TYPE == 'type T':  # true_deltaV13 is different for type T and R
        true_deltaV13 = [-1*(seebeck_SRM3451(TREF_K) - S_Con)*delta_T for
                         delta_T in dT_true]
        true_deltaV24 = [-1*(seebeck_SRM3451(TREF_K) - S_Cu)*delta_T for
                         delta_T in dT_true] #dT_true is a global variable
        # note: true_deltaV is in uV
        # introduce voltage offset for true_deltaV lists
        meas_deltaV13 = [volt + offs_inputs[1] + offs_inputs[3] for volt in true_deltaV13] # uV
        meas_deltaV24 = [volt + offs_inputs[2] + offs_inputs[4] for volt in true_deltaV24] # uV
        # choose between deltaV13 or deltaV24 for Seebeck voltage
        meas_deltaV = [meas_deltaV24, meas_deltaV13] # uV

        use_top_13_wires = False

    elif THERMOCOUPLE_TYPE == 'type R':
        true_deltaV13 = [-1*(seebeck_SRM3452_SiGe(TREF_K) - S_Pt)*delta_T for
                         delta_T in dT_true]
        # note: true_deltaV is in uV
        # introduce voltage offset for true_deltaV lists
        meas_deltaV13 = [volt + offs_inputs[1] + offs_inputs[3] for volt in true_deltaV13] # uV
        # meas_deltaV24 = [volt + offs_inputs[2] + offs_inputs[4] for volt in true_deltaV24] # uV
        meas_deltaV = [None, meas_deltaV13] # uV
        # choose between deltaV13 or deltaV24 for Seebeck voltage
        # NOTE: always use 13 for type R because we don't have S for Pt-13%Rh
        use_top_13_wires = True
    else:
        raise ValueError('global constant THERMOCOUPLE_TYPE is incorrect')

    # get a dictionary with slope, intercept, and trendline y values
    trend_info = calculate_trendline(meas_dT, meas_deltaV[use_top_13_wires])

    # plot DV vs DT after offset affects are calculated
    plt.plot(meas_dT, meas_deltaV[use_top_13_wires], 'k.')
    plt.plot(meas_dT, trend_info['trendline'], 'k', label="including offsets")
    # plt.title('Thermoelectric Votlage Produced by Seebeck Effect in Bi₂Te₃₊ₓ',
    #           pad=20)

    # plot expected DV vs DT with no offsets
    trend_true = calculate_trendline(dT_true, true_deltaV13)

    plt.plot(dT_true, true_deltaV13, 'k.')
    plt.plot(dT_true, trend_true['trendline'], 'k--', label="No offsets")


    plt.title('Thermoelectric Votlage Produced by Seebeck Effect',
              pad=20)
    plt.xlabel('$\Delta$ Temperature (K)', fontsize=font_size)
    plt.ylabel('Seebeck Voltage (uV)', fontsize=font_size)
    plt.legend(bbox_to_anchor=(1.05,1))
    invoke_plot_params("Seebeck_True_Measured_for_Slope")


if enable_examine_blip_region:
    # conditions chosen for closely examining blip:
        # type R thermocouple
        # dV4 = 200uV, Tref = 500K, offset varying in wire 3 from -20 to 20 uV
    offs_inputs = [0, 0, 0, 0, 200]  # units: uV
    # dv4_varying = np.linspace(-20, 20, NUM_POINTS)  # units: uV
    dv4_varying = np.linspace(-5.9, -5.7, NUM_POINTS)  # units: uV
    delta_V34_true = [temp_to_voltage(temp, TREF_C) for temp in Tcolds_C] # mV

    delta_V12_true = [temp_to_voltage(temp, TREF_C) for temp in Thots_C] # mV
    delta_V12_meas = [volt + offs_inputs[1]/1000 + offs_inputs[2]/1000 for volt
                      in delta_V12_true]
    offs_Thots_C = [voltage_to_temp(volt, TREF_C) for volt in delta_V12_meas]

    S_Pt = seebeck_platinum(TREF_K) # units: uV/K
    true_deltaV13 = [-1*(seebeck_SRM3452_SiGe(TREF_K) - S_Pt)*delta_T for
                         delta_T in dT_true]

    dV34_varying = []
    offs_Tcolds_varying = []
    meas_dT_varying = []
    meas_dV13_varying = []
    DV_DT_slope_varying = []

    for voltage in dv4_varying:
        offs_inputs[3] = voltage
        delta_V34_meas = [volt + offs_inputs[3]/1000 + offs_inputs[4]/1000 for
                          volt in delta_V34_true]
        dV34_varying.append(delta_V34_meas)

        offs_Tcolds_C = [voltage_to_temp(volt, TREF_C) for volt in delta_V34_meas]
        offs_Tcolds_varying.append(offs_Tcolds_C)

        meas_dT = [offs_Thots_C[ind]-offs_Tcolds_C[ind] for
                  ind in range(len(offs_Thots_C))]
        meas_dT_varying.append(meas_dT)

        meas_deltaV13 = [volt + offs_inputs[1] + offs_inputs[3] for volt in
                         true_deltaV13] # uV
        meas_dV13_varying.append(meas_deltaV13)

        trend_info = calculate_trendline(meas_dT, meas_deltaV13)
        DV_DT_slope_varying.append(trend_info['slope'])


    plt.plot(dv4_varying, dV34_varying, 'k')
    plt.title('Measured dV34 (mV)', pad=20)
    plt.xlabel('$\delta$V3 ($\mu$V)', fontsize=font_size)
    plt.ylabel('Measured dV34 (mV)', fontsize=font_size)
    invoke_plot_params("meas_dV34")

    plt.plot(dv4_varying, offs_Tcolds_varying, 'k')
    plt.title('Tcolds Including Offsets (C)', pad=20)
    plt.xlabel('$\delta$V3 ($\mu$V)', fontsize=font_size)
    plt.ylabel('Tcolds Including Offsets (C)', fontsize=font_size)
    invoke_plot_params("offs_Tcolds_varying")

    plt.plot(dv4_varying, meas_dT_varying, 'k')
    plt.title('Measured Temperature Difference (C)', pad=20)
    plt.xlabel('$\delta$V3 ($\mu$V)', fontsize=font_size)
    plt.ylabel('Measured Temperature Difference (C)', fontsize=font_size)
    invoke_plot_params("meas_dT_varying")

    plt.plot(dv4_varying, meas_dV13_varying, 'k')
    plt.title('Measured Seebeck Voltage 13 (uV)', pad=20)
    plt.xlabel('$\delta$V3 ($\mu$V)', fontsize=font_size)
    plt.ylabel('Measured Seebeck Voltage 13 (uV)', fontsize=font_size)
    invoke_plot_params("meas_dV13_varying")

    plt.plot(dv4_varying, DV_DT_slope_varying, 'k')
    plt.title('DV vs DT slope (uV/K)', pad=20)
    plt.xlabel('$\delta$V3 ($\mu$V)', fontsize=font_size)
    plt.ylabel('DV vs DT slope (uV/K)', fontsize=font_size)
    invoke_plot_params("DV_DT_slope_varying")


# things that need are changed when switching between type T and type R:
    # thermocouple conversion polynomials
    # seebeck coefficient subtracted out in seebeck measurement function
    # SRM material used (Bismuth Telluride or Silicon 80 Germanium 20)
    # temperature / voltage ranges used to plot polynomials

