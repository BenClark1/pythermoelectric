# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 17:35:54 2021

@author: benja
"""
# import pandas as pd # for plotting
# import matplotlib.pyplot as plt

# import xlsxwriter #can't download xlsxwriter
  
# workbook = xlsxwriter.Workbook('test.xlsx')
# worksheet = workbook.add_worksheet()

# import xlwt
# from xlwt import Workbook
  
# # Workbook is created
# wb = Workbook()
  
# # add_sheet is used to create sheet.
# sheet1 = wb.add_sheet('Sheet 1')

# import pandas as pd

import math

test_x_vals = [-50,
-45,
-40,
-35,
-30,
-25,
-20,
-15,
-10,
-5,
0,
5,
10,
15,
20,
25,
30,
35,
40,
45,
50
]

test_y_vals = [-37.89598731,
-33.87475956,
-28.38121294,
-23.01716199,
-17.11340011,
-12.22167504,
-8.6646922,
-3.53805506,
2.232244935,
7.59703185,
14.09501266,
18.22373359,
22.7946269,
27.5336017,
34.07886067,
39.54899398,
44.60424249,
47.78676679,
51.27480414,
56.90450509,
62.12596717]

data = 0
i = 0
# df_dict = {'Uneg (uV)':x_vals}

def calculate_trendline(x_vals, y_vals):
    # for ind in range(len(x_vals)):
    #     # current_ind = x_vals.index(uvolt)
    #     print(x_vals[ind], end='')
    #     print('\t', end='')
    #     print(y_vals[ind])
        # sheet1.write(current_ind, 4, uvolt)
        # sheet1.write(current_ind, 5, y_vals[current_ind])
    
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
    # correlation /= (sdev_x*sdev_y) #doesn't agree with excel
    # correlation /= (len(x_vals) - 1)
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


# calculate_trendline(test_x_vals, test_y_vals)
# df = pd.DataFrame(df_dict)
# df.to_excel('./test.xlsx', startrow=0, startcol=0)


















