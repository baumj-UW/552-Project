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
from scipy.integrate import odeint  #refs odeint directly instead of long pointer
import matplotlib.pyplot as plt  #refs this pointer as plt --> try simplifiying this later


# simple ODE practice
def model(y,t):
    k1 = 0.3
    k2 = 2.0
    x1, x2 = y
    dydt = [x2, -k1*x2 - k2*numpy.sin(x1)]
    return dydt

#Gen Model --> expand each equation to be a vector representing each gen; or loop function for each gen
def genModel (gens, t):
    swingEqn, deltaW = gens  # vector of variable outputs
    M = 0.6 #2*H*S/w_s 
    Pm = 1.999 
    Pe = 1.75
    DAMP = 1.0 ## does this need to be here?
#     swingEqn = (1/M)*Pm - Pe - DAMP*deltaW # M*delta(w') = Pm - Pe - D*delta(w)
#     speedEqn = deltaW
    dgdt = [(1/M)*Pm - Pe - DAMP*deltaW, deltaW] #, emfEqn]
    return dgdt # dgdt = synch gen model diff eqs

# init condits
gen0 = [0.0, 1.0]

#time points
t = numpy.linspace(0,1.5)  #change time steps 

#solve gen eqns
response = odeint(genModel,gen0,t)
# #plot
plt.plot(t,response[:,0])
plt.plot(t,response[:,1])
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()

# #init condits
# y0 = [numpy.pi -0.1, 0.0]
# 
# #time points
# t = numpy.linspace(0,10,101)
# 
# #solve ode
# sol = odeint(model,y0,t)
# 
# #plot
# plt.plot(t,sol[:,0])
# plt.plot(t,sol[:,1])
# plt.xlabel('time')
# plt.ylabel('y(t)')
# plt.show()


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

# delta(w) = 'speed deviation'
#===============================================================================

 

