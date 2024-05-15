# -*- coding: utf-8 -*-
"""
Created on Mon May 13 10:01:43 2024

@author: heavyhalogen
"""
# Import packages

import numpy as np
import matplotlib.pyplot as plt

# 1.Start
# File parameter check

print('This script applies Raman spectra correction on the water bands and calculates the R\u1d05 value after Grützner and Bureau, 2024.')
print('The dataset has to be a single spectrum from a .txt, .csv, or a related type of flat file.')
print('The shift on the x-axis should be in cm\u207b\u00b9.')

file_name = input('File path and name: ')
rows = int(input('How many rows do you want to skip in your file from the top on?'))
print('Please skip all rows at the top of the file that contain header or comment text.')
limit = input('What is the delimiter? For Tab use backslash+t.')

# Open file and create numpy arrary 
       
with open(file_name, 'r') as file:
    data = np.loadtxt(file, delimiter=limit, dtype=float, skiprows=rows)


# Choosing x and y columns

rd_pos = input('Do you want to set column position of shift and intensity (y/n)? Default is first (1) value = intensity, second (2) value = shift.')
if rd_pos in ('y', 'Y'):
    pos_x = int(input('Intensity: ')) - 1
    pos_y = int(input('Shift: ')) - 1
    x = data[:,pos_x]
    y = data[:,pos_y]
else:
    x = data[:,0]
    y = data[:,1]

    
# Plot the raw data

plt_x_low = 2900
plt_x_high = 3700
x_bar = False

rd_plot = input('Do you want to plot the uncorrected raw spectrum (y/n)?')

if rd_plot in ('y','Y'):
    x_bar = input('Do you want to set the lower and upper limit for the x-axis (y/n)? Default is 2900 to 3700 cm\u207b\u00b9.')
    if x_bar in ('y', 'Y'):
        plt_x_low = float(input('Lower limit of the x-axis: '))
        plt_x_high = float(input('Upper limit of the x_axis: '))
    plt.plot(x,y)
    plt.xlim(plt_x_low, plt_x_high)
    plt.xlabel('cm\u207b\u00b9')
    plt.ylabel('Intensity')
    plt.show





# 2. Linear baseline correction

# 2.1. Find minima at start and end of the range of interest

difference_2890 = np.abs(x - 2890)         
x_bsl_st_2890 = int(difference_2890.argmin())                                    #find lower index value for minimum in start range

difference_3090 = np.abs(x - 3090)         
x_bsl_st_3090 = int(difference_3090.argmin())                                    #find upper index value for minimum in start range

bsl_start_v = min(y[x_bsl_st_2890:x_bsl_st_3090])                                   #Find minimum value's index in end range
bsl_start_t = np.where(y == bsl_start_v)                                        # numpy.where() often gives tuples - but not always.
bsl_start_i = bsl_start_t[0]                                                   # Extract tuple from array
bsl_start_i = np.array(bsl_start_i[-1]).item()                                  #Extract value from tuple



difference_3590 = np.abs(x - 3590)         
x_bsl_end_3590 = int(difference_3590.argmin())                                  #find lower index value for minimum in end range

difference_3710 = np.abs(x - 3710)         
x_bsl_end_3710 = int(difference_3710.argmin())                                  #find upper index value for minimum in end range

bsl_end_v = min(y[x_bsl_end_3590:x_bsl_end_3710])                                   #Find minimum value's index in end range
bsl_end_t = np.where(y == bsl_end_v)
bsl_end_i = bsl_end_t[0]
bsl_end_i = np.array(bsl_end_i[-1]).item()


    
#2.2. Create and apply linear baseline correction

bsl_steps = bsl_end_i - bsl_start_i                                                 # How many steps

bsl_start = y[bsl_start_i]                                                        #Minimum value in end range
bsl_end = y[bsl_end_i]                                                            #Minimum value in end range

bsl_m = (bsl_end_v - bsl_start_v) / bsl_steps                                        # baseline slope
bsln_n = bsl_start - (x[bsl_start_i] * bsl_m)                                         # y-intercept

for i in y:
    bsl = x * bsl_m + bsln_n
    y_bsl = y - bsl


# Plot the baseline corrected data

rd_plot = input('Do you want to plot the baseline-corrected spectrum (y/n)?')

if rd_plot in ('y','Y'):
    if x_bar not in ('y', 'Y'):
        x_bar = input('Do you want to set the lower and upper limit for the x-axis (y/n)? Default is 2900 to 3700 cm\u207b\u00b9.')
        if x_bar in ('y', 'Y'):
            plt_x_low = float(input('Lower limit of the x-axis: '))
            plt_x_high = float(input('Upper limit of the x_axis: '))
    plt.plot(x,y_bsl)
    plt.xlim(plt_x_low, plt_x_high)
    plt.xlabel('cm\u207b\u00b9')
    plt.ylabel('Intensity')
    plt.show





#2.2. Calculate  RD value

turn_v = np.abs(x - 3325)                                                       
turn_i = int(turn_v.argmin())                                                   # Index of turning point
turn_v = y[turn_i]                                                              # turning point or isobestic point

peak_v = max(y[turn_i:x_bsl_end_3590])                                          # right peak point
peak_t = np.where(y == peak_v)
peak_i = peak_t[0]
peak_i = np.array(peak_i[-1]).item()                                          # Index of right peak point

steps_right = peak_i - turn_i                                                   # steps between turning point and peak


rd = np.sum(y[turn_i:peak_i]) / (turn_v * steps_right)                          # RD value







#2.3. normalization is not applied for single spectra



#2.4. Smoothing

j = 0
y_av = []                                                                       
for i in y_bsl:
    k = j - 20
    l = j + 20
    while k < 0:
        k = 0 
    while l > len(y_bsl) - 1:
        l = len(y_bsl) -1
    av_j = np.average(y_bsl[k:l])                                                   # adjacent averaging with step size 40 (+/- 20)
    y_av.append(av_j)
    j = j + 1

y_sm = np.array(y_av) 
  
# Plot the data

rd_plot = input('Do you want to plot the smoothed spectrum (y/n)?')

if rd_plot in ('y','Y'):
    if x_bar not in ('y', 'Y'):
        x_bar = input('Do you want to set the lower and upper limit for the x-axis (y/n)? Default is 2900 to 3700 cm\u207b\u00b9.')
        if x_bar in ('y', 'Y'):
            plt_x_low = float(input('Lower limit of the x-axis: '))
            plt_x_high = float(input('Upper limit of the x_axis: '))
    plt.plot(x,y_sm)
    plt.xlim(plt_x_low, plt_x_high)
    plt.xlabel('cm\u207b\u00b9')
    plt.ylabel('Intensity')
    plt.show
    print('Note that smoothing is applied after the R\u1d05 calculation.')


# Print RD value at the end

print('')
print('')
print('R\u1d05 value after Grützner and Bureau, 2024: ', np.around(rd, 3))

