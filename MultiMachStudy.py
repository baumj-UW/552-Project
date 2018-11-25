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
from numpy.linalg import inv
from scipy.integrate import odeint  #refs odeint directly instead of long pointer
import matplotlib.pyplot as plt  #refs this pointer as plt --> try simplifiying this later


#Define Constants

NGEN = 3 
NBUS = 8


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
# plt.plot(t,response[:,0])
# plt.plot(t,response[:,1])
# plt.xlabel('time')
# plt.ylabel('y(t)')
# plt.show()

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

# Step 1 - Convert to common base (already done for this project) 

# Create Ybus matrix
 
Ybus_pre = numpy.array([[-12.5j, 0.0, 0.0, 12.5j, 0.0, 0.0, 0.0, 0.0],
                    [0.0, -5.556j, 0.0, 0.0, 5.556j, 0.0, 0.0, 0.0],
                    [0.0, 0.0, -8.333j, 0.0, 0.0, 8.333j, 0.0, 0.0],
                    [12.5j, 0.0, 0.0, -32.48j, 10.0j, 0.0, 10.0j, 0.0],
                    [0.0, 5.556j, 0.0, 10.0j, -35.526j, 0.0, 10.0j, 10.0j],
                    [0.0, 0.0, 8.333j, 0.0, 0.0, -28.0j, 10.0j, 10.0j],
                    [0.0, 0.0, 0.0, 10.j, 10.0j, 10.0j, 2.917-31.217j, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 10.0j, 10.0j, 0.0, 1.363-20.369j]], dtype = complex)  #append definition with dtype=complex

Ybus_fault = numpy.array([[-12.5j, 0.0, 0.0, 12.5j, 0.0, 0.0, 0.0, 0.0],
                    [0.0, -5.556j, 0.0, 0.0, 5.556j, 0.0, 0.0, 0.0],
                    [0.0, 0.0, -8.333j, 0.0, 0.0, 8.333j, 0.0, 0.0],
                    [12.5j, 0.0, 0.0, -32.48j, 10.0j, 0.0, 0.0, 0.0],
                    [0.0, 5.556j, 0.0, 10.0j, -35.526j, 0.0, 0.0, 10.0j],
                    [0.0, 0.0, 8.333j, 0.0, 0.0, -28.313j, 0.0, 10.0j],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 10.0j, 10.0j, 0.0, 1.363-20.369j]], dtype = complex) 

Ybus_post = numpy.array([[-12.5j, 0.0, 0.0, 12.5j, 0.0, 0.0, 0.0, 0.0],
                    [0.0, -5.556j, 0.0, 0.0, 5.556j, 0.0, 0.0, 0.0],
                    [0.0, 0.0, -8.333j, 0.0, 0.0, 8.333j, 0.0, 0.0],
                    [12.5j, 0.0, 0.0, -32.48j, 10.0j, 0.0, 10.0j, 0.0],
                    [0.0, 5.556j, 0.0, 10.0j, -35.526j, 0.0, 10.0j, 10.0j],
                    [0.0, 0.0, 8.333j, 0.0, 0.0, -18.313j, 0.0, 10.0j],
                    [0.0, 0.0, 0.0, 10.j, 10.0j, 0.0, 2.917-21.217j, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 10.0j, 10.0j, 0.0, 1.363-20.369j]], dtype = complex) 

# Step 2 Add model of load admittances

# Step 3 calculate internal gen voltages

# Step 4 cacluate prefault / fault / post fault admittance matrices 

# Step 5 Kron reduction 
# Yreduced = Ynn - Yns*(1/Yss)*Ysn
# n = NGENS , s = remaining system nodes
def kronRed(Y,n,s): 
    Ynn = Y[0:n,0:n]
    Yns = Y[0:n,n:s] 
    Ysn = Y[n:s,0:n]
    Yss_inv = inv(Y[n:s,n:s])
    Yhat = Ynn - numpy.dot(numpy.dot(Yns,Yss_inv),Ysn)
    return Yhat 

Ypre_red = kronRed(Ybus_pre,NGEN,NBUS)
#Yfault_red = kronRed(Ybus_fault,NGEN,NBUS)  --> fix singular issue with fault bus
Ypost_red = kronRed(Ybus_post,NGEN,NBUS)
print(Ypre_red)
#print(Yfault_red)
print(Ypost_red)
#===============================================================================
# Solve Equations with scipy.integrate.odeint 
# M*delta(w') = Pm - Pe - D*delta(w)
# d' = delta(w)
# Ttransdo*Etransq' = Ef - Etransq + Id*(Xd - Xtransd)

# delta(w) = 'speed deviation'
#===============================================================================

 

