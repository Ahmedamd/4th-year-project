# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 09:45:30 2018

@author: amand and ahmed
"""

import stamps
import copy
import numpy as np
import math
from numpy.linalg import inv
import sympy

#netlist = str('freq,10,50,80,V,0 -9,GND,101,R,0 -9,5 -9,0.1,R1,0 -9,7 -20,50,L,5 -9,GND,0.000001,C,5 -9,GND,100')
#netlist = ('freq,10,10,500, V,-11 -4,GND,5,R,-11 -4,-6 -4,1000,C,-6 -4,-1 -4,0.001,L,GND,-1 -4,0.00001,R,GND,-6 -4,1000')
#netlist = ('freq,1,10000,10000, V,-7 -7,GND,5,R,-7 -7,-2 -7,1000,R,-2 -7,3 -7,1000,R,GND,-2 -7,1000,C,3 -7,GND,0.001')
netlist = ('time,1,100,5, sine,-7 -7,GND,100000,R,-7 -7,-2 -7,100,sine,-2 -7,GND, 500')
print(netlist)

words = netlist.split(",") 
print(words)
k = []
store = []
count = 1
for element in words:
    store.append(element)
    if count == 4:
        count = 1
        k.append(store)
        store = []
    else:
        count = count + 1
        
print("kkkkkkk\n", k)
        
        
#import os
#import re
#specialchar = ';'
#k = ['R,-13 -4,-8 -4,1000','L,-8 -4,-2 -4,0.00001','L,-2 -4,4 -4,0.00001']
##k = k.partition(specialchar)
#print (k)
##k = k.splitlines()
dict = {}#creating a dictionary
new_matrix = [] #manipulated values will be stored in matrix
counter = 0
current_top = 0
temp = 0
i=0
for line in k:
    if i == 0:
        new_matrix.append(line)
        i=i+1
    else:
        print(line)
        internal_arr = []
        comp_start = line[0]
        node1 = line[1]
        node2 = line[2]
        comp_end = line[3]
        print("com_end", comp_end)
        internal_arr.append(comp_start)
        if node1 not in dict:
            if counter == 0 and node1 == 'GND':
                dict[node1] = 0
            elif node1 == 'GND':
                dict[node1] = 0
     
            else:
                temp = counter + 1
                dict[node1] = temp
        if node2 not in dict:
            if node2 == 'GND':
                dict[node2] = 0
            else:
                dict[node2] = counter + 1
        internal_arr.append(dict[node1])
        internal_arr.append(dict[node2])
        internal_arr.append(comp_end.strip())
        new_matrix.append(internal_arr)
        counter += 1
        
        print("internal \n", internal_arr)
print (new_matrix)

# Entering frequency domain mode
if 'freq' in new_matrix[0][0] : 
    matrixA = []
    matrixG = []
    matrixC = []
    matrixX = []
    matrixZ = []
    resistors = []
    capacitors = []
    inductors = []
    voltage_source = []
    current_source = []
    circuit_element = []
    point1 = []
    freq_range = []
    
    #m & n will allow us to find the dimensions of the matrix
    n = 0
    m = 0
    b = 0
    
    index = 0
    
    
    for element in new_matrix:
     
        print(element)
        
    
        if index > 0 :
            #Finding out the number of nodes
            if (int(element[1]) > n):
                n = int(element[1])
            if (int(element[2]) > n):
                n = int(element[2])
                
    
        index = index + 1
        #Placing all circuit elements in the appropriate array
        if "R" in element[0]:
            resistors.append(element)
        if "C" in element[0]:
            capacitors.append(element)
        if "L" in element[0]:
            inductors.append(element)
            b=b+1
        if "I" in element[0]:
            current_source.append(element)
        if "V" in element[0]:
            voltage_source.append(element)
            m=m+1
        if "f" in element[0]:
            freq_range.append(element)
            
    
    
    
        #creating an empty matrix A with the correct dimensions
        matrixG = np.zeros((m+n+b,m+n+b),dtype=np.complex)
        #creating an empty matrix C for the reactive components with the correct dimensions
        matrixC = np.zeros((m+n+b,m+n+b),dtype=np.complex)
        #creating an empty matrix Z with the correct dimensions
        matrixZ = np.zeros((m+n+b,1),dtype=np.complex)
    
    
    plots = []
    phase_plots = []
    
    frequency = stamps.F(freq_range)
    
    
    length = len(frequency)
    
    print(length)
    plots = np.zeros((m+n+b,length))
    phase_plots = np.zeros((m+n+b,length))
    stamps.R(resistors, matrixG)
    stamps.V(voltage_source, matrixG, matrixZ, n, b)
    stamps.I(current_source, matrixZ)
    
    
    stamps.C(capacitors, matrixC)
    stamps.L(inductors, matrixG, matrixC, n)
    
    #print("\n matrixG", matrixG, "\n")
    #print("\n matrixC", matrixC, "\n")
    
    index = 0
    
    print("\n STARRRRTTTTT \n\n\n")
    for element in frequency:
        
        s=complex(1j*2*pi*element)
        
      
        #creating an empty matrix A with the correct dimensions
        #We want this matrix to reset to 0 at the beginning of each loop
        matrixA = np.zeros((m+n+b,m+n+b),dtype=np.complex)
        #creating an empty matrix X with the correct dimensions
        matrixX = np.zeros((m+n+b,1), dtype=np.complex)
        
        #print(s)
        #print(s*matrixC)
        
        
        matrixA = matrixG + (s*matrixC)
        
        
    #    print("ELEMENT: ", element)
    #    print("\n matrixG:\n", matrixG, "\n" )
    #    print("\n matrixC:\n", matrixC, "\n" )
    #    print("\n matrixA:\n", matrixA, "\n" )
    #    print("\n matrixZ:\n", matrixZ, "\n" )
    #    
    #    print("matrixA\n", matrixA)
    #    print("invA\n", invA)
    #    
        
        newA = matrixA
        
    
        matrixX = np.linalg.solve(matrixA, matrixZ)
    
       # print("\n matrixA\n", plots, "\n")
        #print(matrixX)
    #    print("\n matrixX:\n", matrixX)
        
        stamps.Mag_points(plots, matrixX, index, m, n, b)
        stamps.Phase_points(phase_plots, matrixX, index, m, n, b)
        
    #    plots.append(magnitude1)  
        
        #print(matrixX)
        #print("\n", abs(matrixX[3]), "\n")
        index = index+1
    
    #    print(plots)
    #print(matrixA)
    #
    for i in range(0, (n+m)):
        
        
    #    print("radians attempt")
    #    plt.plot(2*180*frequency, plots[i], color='blue')
    #    plt.xscale('log')
    #    plt.xlabel("Frequency (Hz)")
    #    plt.ylabel("Magnitude (dB)")
    #    plt.title("Node %i" %(i+1))
    #    plt.grid()
    #    plt.show()
        
        print("Magnitude")
        plt.subplot(n+m,1,i+1)
        plt.plot(frequency, plots[i], color='blue')
        plt.xscale('log')
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Magnitude (dB)")
        plt.title("Node %i" %(i+1))
        plt.grid()
        plt.show()
        
        print("Phase")
        plt.subplot(n+m,1,i+1)
        plt.plot(frequency, phase_plots[i], color='blue')
        plt.xscale('log')
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Degrees")
        plt.title("Node %i" %(i+1))
        plt.grid()
        plt.show()

    
    
# Entering time domain mode
if 'time' in new_matrix[0][0]:
    matrixA = []
    matrixG = []
    matrixC = []
    matrixX = []
    matrixZ = []
    resistors = []
    capacitors = []
    inductors = []
    voltage_source = []
    current_source = []
    sine = []
    time_range = []
    time_plots = []
    sine_steps = []
    time = []
    
    #m & n will allow us to find the dimensions of the matrix
    n = 0
    m = 0
    b = 0
    s = 0
    
    
    
    #Reading the netlist lines and placing them in an array.
    #Each array element now corresponds to a circuit element
    #And the information corresponding to it (i.e nodes)
    
    
    
    
    
    index = 0
    for element in new_matrix:
        
        
        if index > 0 :
        
            #Finding out the number of nodes
            if (int(element[1]) > n):
                n = int(element[1])
            if (int(element[2]) > n):
                n = int(element[2])
                    
            
        index = index + 1
        #Placing all circuit elements in the appropriate array
        if "R" in element[0][0]:
            resistors.append(element)
        if "C" in element[0][0]:
            capacitors.append(element)
        if "L" in element[0][0]:
            inductors.append(element)
            b=b+1
        if "I" in element[0][0]:
            current_source.append(element)
        if "V" in element[0][0]:
            voltage_source.append(element)
            m=m+1
        if "t" in element[0][0]:
            time_range.append(element)
        if "s" in element[0][0]:
            sine.append(element)
            s = s+1
    
    
    
    #Obtaining that time range of the analysis
    #start = Decimal(time_range[0][1].strip('"'))
    #stop = Decimal(time_range[0][2].strip('"'))
    #step_size = Decimal(time_range[0][3].strip('"')) 
    #samples = math.floor((stop - start)/step_size) + 1
    
    start = float(time_range[0][1])
    stop = float(time_range[0][2])
    step_size = float(time_range[0][3]) 
    samples = math.floor((stop - start)/step_size) + 1
    
    print(samples)
    
    
    
    #creating an empty matrix G with the correct dimensions 
    #the G matrix contains resistors and the 1/-1 of voltage sources or inductors
    matrixG = np.zeros((m+n+b+s,m+n+b+s))
    
    #creating an empty matrix C for the reactive components with the correct dimensions
    matrixC = np.zeros((m+n+b+s,m+n+b+s))
    
    #creating an empty matrix z with the correct dimensions
    matrixZ = np.zeros((m+n+b+s,1))
    
    #creating an empty matrix with the correct dimensions
    matrixRightSide = np.zeros((m+n+b+s,1))
    
    #creating an empty matrix A with the correct dimensions
    matrixA = np.zeros((m+n+b+s,1))
    
    # matrix of the unknowns
    matrixX = np.zeros((m+n+b+s, 1))
    
    newmatrixX = np.zeros((m+n+b+s, 1))
    
    #This matrix will contain all of the data points that need to be plotted
    # ???????????????????????? ask Pavan about(m+n+b+s)
    time_plots = np.zeros((m+n+b+s, samples))
    
    # array containing all the time instances to be plotted
    
    #Creating a matrix to contain the time increments/steps/values for the sine generators
    t =np.zeros((1, int(samples)))
    t = t[0]
    
    #creating an array cointaining all of the points from the sine
    #wave generator corresponding to the desired steps
    sine_steps = np.zeros((s, samples))
    
    
    #Inserting the resistors, voltages(not sine or pule or etc..), current,
    #capacitors, inductors
    stamps.R_TD(resistors, matrixG)
    stamps.V_TD(voltage_source, sine, matrixG, matrixZ, n, b, m)
    stamps.I_TD(current_source, matrixZ)
    stamps.C_TD(capacitors, matrixC)
    stamps.L_TD(inductors, matrixG, matrixC, n)
    
    stamps.sine_generator_TD(start, stop, step_size, sine_steps, sine, samples, s, t)
    
    
    # Initiating the time domain calculations
    # First step is a special case: 
    
    
    index = 0
    for i in range(s):
        matrixZ[m+n+b+i] = sine_steps[i][0]
        print("initial", sine_steps[i][0])
        matrixX = np.linalg.solve(matrixG, matrixZ)
        stamps.Time_points_TD(time_plots, matrixX, index, m, n, b, s) 
        
        print(matrixG)
        print(matrixZ)
        print(matrixX)
        print(time_plots)
        
    # Calculating for all of the other instances
        
    # calculating matrix A = (G + C/h)
    matrixA =  matrixG + np.divide(matrixC, float(time_range[0][3]))
    
    
    index = 1
    
    for j in range(samples-1):
        
        for i in range(s):
            matrixZ[m+n+b+i] = sine_steps[i][index]
            matrixRightSide = matrixZ + np.divide(matrixC, float(time_range[0][3])) * matrixX
            matrixX = np.linalg.solve(matrixA, matrixRightSide)
            
            
            for k in range(m+n+b+s):
                #THIS IS WHERE WE LEFT OFF
                #MATRIX X CURRENTLY HAS THE WRONG DIMENSIONS> 
                #NEEDS TO BE FIXED
                newmatrixX[k] = matrixX[k][0]
                
            stamps.Time_points_TD(time_plots, newmatrixX, index, m, n, b, s) 
            print("matrixX \n", matrixX)
            print("newmatrixX \n", newmatrixX)
            
        index = index + 1    
    
    for i in range(0, (n+m+b+s)):
        
        
        print("radians attempt")
        plt.plot(t, time_plots[i], color='blue')
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Magnitude (dB)")
        plt.title("Node %i" %(i+1))
        plt.grid()
        plt.show()