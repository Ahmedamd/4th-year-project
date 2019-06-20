# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 22:06:10 2018

@author: amand and ahmed
"""

import numpy as np
from numpy import pi
import math
import matplotlib.pyplot as plt

def R(resistors, matrixG):
    #Resistor stamps  
    node1=0
    node2=0
    for element in resistors:
        
        
        node1 = int(element[1])
        node2 = int(element[2])
        if node1>0:
            matrixG[node1-1, node1-1] = matrixG[node1-1, node1-1] + 1/int(element[3])
        if node2>0:
            matrixG[node2-1, node2-1] = matrixG[node2-1, node2-1] + 1/int(element[3])
            if node1 >0:
                matrixG[node1-1, node2-1] = matrixG[node1-1, node2-1] - 1/int(element[3])
                matrixG[node2-1, node1-1] = matrixG[node2-1, node1-1] - 1/int(element[3])
                
        #print(matrixA)
    return(matrixG)
  
    

def V(voltage_source, sine, matrixG, matrixZ, n, b, m):
    #Voltage sources        
    #remember that node1 is always connected to the
    #positive terminal and node2 to the negative terminal
                                                    
    index = 0
    for element in voltage_source:
        node1 = int(element[1])
        node2 = int(element[2])
        if node1>0:
            matrixG[node1-1, n+b+index] = float(1)
            matrixG[n+b+index, node1-1] = float(1)
        if node2>0:
            matrixG[n+b+index, node2-1] = float(-1)
            matrixG[node2-1, n+b+index] = float(-1)
        index = index+1
        
    count = 0
    for element in voltage_source:
        matrixZ[n+b+count] = matrixZ[n+b+count] + float(element[3])
        count = count+1
        
    index = 0
    for element in sine:
        print("sffsf")
        node1 = int(element[1])
        node2 = int(element[2])
        if node1>0:
            matrixG[node1-1, n+b+m+index] = float(1)
            matrixG[n+b+m+index, node1-1] = float(1)
        if node2>0:
            matrixG[n+b+m+index, node2-1] = float(-1)
            matrixG[node2-1, n+b+m+index] = float(-1)
        index = index+1
        
    return(matrixG, matrixZ)




def C(capacitors, matrixC):
    #Calculate matrix for every frequency. 
   
    node1=0
    node2=0
    
    for element in capacitors:
        node1 = int(element[1])
        node2 = int(element[2])

        if node1>0:
            matrixC[node1-1, node1-1] = matrixC[node1-1, node1-1] + float(element[3])
        if node2>0:
            
            matrixC[node2-1, node2-1] = matrixC[node2-1, node2-1] + float(element[3])
            if node1 >0:
                matrixC[node1-1, node2-1] = matrixC[node1-1, node2-1] - float(element[3])
                matrixC[node2-1, node1-1] = matrixC[node2-1, node1-1] - float(element[3])
                
    return(matrixC)


def L(inductors, matrixG, matrixC, n):
    #Calculate matrix for every frequency.
    
    node1=0
    node2=0
    index = 0
   
    for element in inductors:
        node1 = int(element[1])
        node2 = int(element[2])
        
        if node1>0:
            matrixG[node1-1, n+index] = float(1)
            matrixG[n+index, node1-1] = float(1)
            
        if  node2>0:
            matrixG[node2-1, n+index] = float(-1)
            matrixG[n+index, node2-1] = float(-1)
        
        matrixC[n+index, n+index] = -(float(element[3]))
        
        
        index = index+1                
                
    return(matrixC, matrixG)


def I(current_source, matrixZ):
    
    node1=0
    node2=0
    for element in current_source:
        node1 = int(element[1])
        node2 = int(element[2])
        
    if node1>0:
        matrixZ[node1-1] = matrixZ[node1-1] + int(element[3])
    if node2>0:
        matrixZ[node2-1] = matrixZ[node2-1] - int(element[3])
        
    return(matrixZ)

def B(matrixB, sine_steps, m, n, b, s):    
    for i in range(s):
        matrixB[m+n+b+i] = sine_steps[i][0]
        print(sine_steps[i][0])
        
        

def sine_generator(start, stop, step_size, sine_steps, sine, samples, s, t):
    

    
    index = 0
    samples = float(samples)
    for element in sine:
        
        f = int(element[3])
        w = 2 * pi * f
        
        
        for i in range(int(samples)):
          
            if i == 0:
                t[i] = start
            else:
                t[i] = t[i-1] + step_size
                
        y = np.sin(w * t)
        
        k = 0    
        for element in y:
            
            sine_steps[index][k] = y[k]
            #print("sine steps\n", sine_steps)
            
            k = k+1
        

        index = index + 1
        
        #print("y", y)
#    print("sine steps\n", sine_steps)
#    print("t", t)
    
#    for i in range(int(s)):
#        
#        plt.plot(t, sine_steps[i], color='blue')
#        plt.grid()
#        plt.show()        
            
        
        
def Time_points(time_plots, matrixX, index, m, n, b, s):
    
    for i in range((m+n+b+s)):
        time_plots[i][index] = matrixX[i] 
        
        
        
        
    
        


        
        
       