'''
Created on Nov 18, 2018

@author: baumj

Multimachine Stability Study 
Power System Dynamics - EE 552

Based on Bergen and Vittal Example 14.6, p. 570
'''

import math     # Math functions 
import cmath    # Complex math function conj, rect
import openpyxl # Methods to read and write xlsx files
import numpy    # Methods for linear algebra

# Create Ybus matrix

Ybus = numpy.array([[-12.5, 0, 0, 12.5, 0, 0, 0, 0],
                    [0, -5.556, 0, 0, 5.556, 0, 0, 0],
                    [0, 0, -8.333, 0, 0, 8.333, 0, 0],
                    [12.5, 0, 0, -32.48, 10, 0, 10, 0],
                    [0,]
                    [0,]
                    [0,]
                    [0,]])  #append definition with dtype=complex


#someone came in and edited the master file!
#panic, how do these get merged?

#make more changes
