# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 22:41:25 2022

@author: benja
"""

# Define coefficients for the thermocouple polynomial conversions




# coefficients for EMF as a function of temperature b0-b14 (-270C<=T<=0C)
b_neg_typeT = [0.0,
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
b_pos_typeT = [0.0,
0.387481063640e-1,
0.332922278800e-4,
0.206182434040e-6,
-0.218822568460e-8,
0.109968809280e-10,
-0.308157587720e-13,
0.454791352900e-16,
-0.275129016730e-19]
# b = [b_neg, b_pos] # b[1] gives positive polynomial, b[0] gives negative
# coefficients for temperature as a function of EMF c0-7 (-200C<=T<=0C)
c_neg_typeT = [0.0,
2.5949192e1,
-2.1316967e-1,
7.9018692e-1,
4.2527777e-1,
1.3304473e-1,
2.0241446e-2,
1.2668171e-3] 
# coefficients for voltage to temp for positive temp/voltage (0C<=T<=400C)
c_pos_typeT = [0.0,
2.592800e1,
-7.602961e-1,
4.637791e-2,
-2.165394e-3,
6.048144e-5,
-7.293422e-7]
# c = [c_neg, c_pos] # c[1] gives positive polynomial, c[0] gives negative

# --------------------------------------------------------

# define type R thermocouple polynomial coefficients (from NIST) 
b_r1_typeR = [0.0, # range 1 coefficients -50C to 1064.18C
5.28961729765e-03, # b represents temp to voltage conversion
1.39166589782e-05,
-2.38855693017e-08,
3.56916001063e-11,
-4.62347666298e-14,
5.00777441034e-17,
-3.73105886191e-20,
1.57716482367e-23,
-2.81038625251e-27]

b_r2_typeR = [2.95157925316e+00, # range 2 coefficients 1064.18C to 1664.5C
-2.52061251332e-03,
1.59564501865e-05,
-7.64085947576e-09,
2.05305291024e-12,
-2.93359668173e-16]

b_r3_typeR = [1.52232118209e+02, # range 3 coefficients 1664.5C to 1768.1C
-2.68819888545e-01,
1.71280280471e-04,
-3.45895706453e-08,
-9.34633971046e-15]

c_r1_typeR = [0.0, # range 1 coefficients -0.226 to 1.923 mV
1.8891380e2, # c represents "reverse" voltage to temp conversion
-9.3835290e1,
1.3068619e2,
-2.2703580e2,
3.5145659e2,
-3.8953900e2,
2.8239471e2,
-1.2607281e2,
3.1353611e1,
-3.3187769e0] 

c_r2_typeR = [1.334584505e1, # range 2 coefficients 1.923 to 13.228 mV
1.472644573e2,
-1.844024844e1,
4.031129726e0,
-6.249428360e-1,
6.468412046e-2,
-4.458750426e-3,
1.994710149e-4,
-5.313401790e-6,
6.481976217e-8]

c_r3_typeR = [-8.199599416e1, # range 3 coefficients 11.361 to 19.739 mV
1.553962042e2,
-8.342197663e0,
4.279433549e-1,
-1.191577910e-2,
1.492290091e-4]

c_r4_typeR = [3.406177836e4, # range 4 coefficients 19.739 to 21.103 mV
-7.023729171e3,
5.582903813e2,
-1.952394635e1,
2.560740231e-1]
