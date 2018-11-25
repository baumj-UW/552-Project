'''
Created on Nov 18, 2018

@author: baumj

Multimachine Stability Study 
Power System Dynamics - EE 552

Based on Bergen and Vittal Example 14.6, p. 570
Machowski pdf p. 486, 3rd order model
'''

import math     # Math functions 
import cmath    # Complex math function conj, rect
import openpyxl # Methods to read and write xlsx files
import numpy    # Methods for linear algebra
from scipy.integrate import odeint


# Create Ybus matrix
# 
# Ybus = numpy.array([[-12.5, 0, 0, 12.5, 0, 0, 0, 0],
#                     [0, -5.556, 0, 0, 5.556, 0, 0, 0],
#                     [0, 0, -8.333, 0, 0, 8.333, 0, 0],
#                     [12.5, 0, 0, -32.48, 10, 0, 10, 0],
#                     [0,]
#                     [0,]
#                     [0,]
#                     [0,]])  #append definition with dtype=complex


# Step 2 Add model of load admittances

# Step 3 calculate internal gen voltages

# Step 4 cacluate prefault / fault / post fault admittance matrices 

# Step 5 Kron reduction 
#===============================================================================
# Solve Equations with scipy.integrate.odeint 
# M*delta(w') = Pm - Pe - D*delta(w)
# d' = delta(w)
# Ttransdo*Etransq' = Ef - Etransq + Id*(Xd - Xtransd)
#===============================================================================

 

