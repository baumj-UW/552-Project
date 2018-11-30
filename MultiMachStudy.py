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
import numpy as np  # Methods for linear algebra
from numpy.linalg import inv
from scipy.integrate import odeint  #refs odeint directly instead of long pointer
from scipy.integrate import solve_ivp #ODE45 equiv (use this instead of odeint)
import matplotlib.pyplot as plt  #refs this pointer as plt --> try simplifiying this later

#Define Constants

NGEN = 3 
NBUS = 8

BUS_CONN = np.array([[1,0,0,1,0,0,0,0],
                    [0,1,0,0,1,0,0,0],
                    [0,0,1,0,0,1,0,0],
                    [1,0,0,1,1,0,1,0],
                    [0,1,0,1,1,0,1,1],
                    [0,0,1,0,0,1,1,1],
                    [0,0,0,1,1,1,1,0],
                    [0,0,0,0,1,1,0,1]]) #T/F matrix indiciating bus connections for the system 

LINE_XL = 0.1 #Line series admittance
LINE_YB = 0.01 #Line shunt admittance 

# simple ODE practice
def model(y,t):
    k1 = 0.3
    k2 = 2.0
    x1, x2 = y
    dydt = [x2, -k1*x2 - k2*np.sin(x1)]
    return dydt

# #Gen Model --> expand each equation to be a vector representing each gen; or loop function for each gen
# def genModel (gens, t):
#     swingEqn, deltaW = gens  # vector of variable outputs
#     M = 0.6 #2*H*S/w_s 
#     Pm = 1.999 
#     Pe = 1.75
#     DAMP = 1.0 ## does this need to be here?
# #     swingEqn = (1/M)*Pm - Pe - DAMP*deltaW # M*delta(w') = Pm - Pe - D*delta(w)
# #     speedEqn = deltaW
#     dgdt = [(1/M)*Pm - Pe - DAMP*deltaW, deltaW] #, emfEqn]
#     return dgdt # dgdt = synch gen model diff eqs
# 
# # init condits
# gen0 = [0.0, 1.0]  #init condits [delta0, w0]
# 
# #time points
# t = np.linspace(0,1.5)  #change time steps 
# 
# #solve gen eqns
# response = odeint(genModel,gen0,t)
# #plot
# plt.plot(t,response[:,0])
# plt.plot(t,response[:,1])
# plt.xlabel('time')
# plt.ylabel('y(t)')
# plt.show()

# #init condits
# y0 = [np.pi -0.1, 0.0]
# 
# #time points
# t = np.linspace(0,10,101)
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


Vmag = np.zeros((NBUS,1))   #vector of voltage magnitudes
Vtheta = np.zeros((NBUS,1))  # vector of volt angles in degrees
Pbus = np.zeros((NBUS,1)) #vector of P at each bus
Qbus = np.zeros((NBUS,1)) #vector of Q at each bus

# cmath.rect(r, phi) --> to convert polar to rect. value to combine Vmag and Vtheta

Xd_trans = np.array([[0.08],[0.18],[0.12]]) #transient reactance from book table

#assign outputs from book solution 
Vmag[NGEN:NBUS] = [[1.04],[1.02],[1.05],[0.9911],[1.0135]]
## Degrees Vtheta[NGEN:NBUS] = [[0.0],[-3.55],[-2.90],[-7.48],[-7.05]]
Vtheta[NGEN:NBUS] = [[0.0],[-0.06196],[-0.05061],[-0.13055],[-0.12305]]

Pbus[0:NGEN] = [[1.9991],[0.6661],[1.600]]
Qbus[0:NGEN] = [[0.8134],[0.2049],[1.051]]
# Create Ybus matrix (without load and gen admittance)

#Calculate equiv. shunt load admittance and add to Ybus matrix
#Create separate G and B buses for real and imag parts
Ybus_pre_G = np.zeros((NBUS,NBUS))
Ybus_pre_B = np.zeros((NBUS,NBUS))

#gen bus series and self admittances
for bus_i in range(0,NGEN):
    for bus_j in range(0,NBUS):
        if BUS_CONN[bus_i,bus_j]:
           if bus_i == bus_j:
               Ybus_pre_B[bus_i,bus_j] = -1/Xd_trans[bus_i]
               Ybus_pre_B[bus_j,bus_i] = -1/Xd_trans[bus_i]
           else:
               Ybus_pre_B[bus_i,bus_j] = 1/Xd_trans[bus_i]
               Ybus_pre_B[bus_j,bus_i] = 1/Xd_trans[bus_i]

#System bus series admittances
for bus_i in range(NGEN,NBUS):
    for bus_j in range(NGEN,NBUS):
        if np.logical_and(BUS_CONN[bus_i,bus_j],(bus_i != bus_j)):
            Ybus_pre_B[bus_i,bus_j]= 1/LINE_XL

#System bus self admittances             
for bus_i in range(NGEN,NBUS):
    Ybus_pre_B[bus_i,bus_i] = -1*np.sum(Ybus_pre_B[bus_i,:]) + LINE_YB*(np.sum(BUS_CONN[bus_i,NGEN:NBUS])-1)



# use complexY = Ybus_pre_G + 1j*Ybus_pre_B
            
#Calculate equiv gen bus E and delta
#Create full augmented Ybus matrix (includes gen and load admittance)
# Pre-fault: calc. from steady state augmented matrix
# fault: Assume 3ph-g fault, set 


#replace this with calculation from actual bus data 
 
Ybus_pre = np.array([[-12.5j, 0.0, 0.0, 12.5j, 0.0, 0.0, 0.0, 0.0],
                    [0.0, -5.556j, 0.0, 0.0, 5.556j, 0.0, 0.0, 0.0],
                    [0.0, 0.0, -8.333j, 0.0, 0.0, 8.333j, 0.0, 0.0],
                    [12.5j, 0.0, 0.0, -32.48j, 10.0j, 0.0, 10.0j, 0.0],
                    [0.0, 5.556j, 0.0, 10.0j, -35.526j, 0.0, 10.0j, 10.0j],
                    [0.0, 0.0, 8.333j, 0.0, 0.0, -28.0j, 10.0j, 10.0j],
                    [0.0, 0.0, 0.0, 10.j, 10.0j, 10.0j, 2.917-31.217j, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 10.0j, 10.0j, 0.0, 1.363-20.369j]], dtype = complex)  #append definition with dtype=complex

#"remove" bus 7 during fault 
Ybus_fault = np.array([[-12.5j, 0.0, 0.0, 12.5j, 0.0, 0.0, 0.0],
                    [0.0, -5.556j, 0.0, 0.0, 5.556j, 0.0, 0.0],
                    [0.0, 0.0, -8.333j, 0.0, 0.0, 8.333j, 0.0],
                    [12.5j, 0.0, 0.0, -32.48j, 10.0j, 0.0, 0.0],
                    [0.0, 5.556j, 0.0, 10.0j, -35.526j, 0.0, 10.0j],
                    [0.0, 0.0, 8.333j, 0.0, 0.0, -28.313j, 10.0j],
                    [0.0, 0.0, 0.0, 0.0, 10.0j, 10.0j, 1.363-20.369j]], dtype = complex) 

Ybus_post = np.array([[-12.5j, 0.0, 0.0, 12.5j, 0.0, 0.0, 0.0, 0.0],
                    [0.0, -5.556j, 0.0, 0.0, 5.556j, 0.0, 0.0, 0.0],
                    [0.0, 0.0, -8.333j, 0.0, 0.0, 8.333j, 0.0, 0.0],
                    [12.5j, 0.0, 0.0, -32.48j, 10.0j, 0.0, 10.0j, 0.0],
                    [0.0, 5.556j, 0.0, 10.0j, -35.526j, 0.0, 10.0j, 10.0j],
                    [0.0, 0.0, 8.333j, 0.0, 0.0, -18.313j, 0.0, 10.0j],
                    [0.0, 0.0, 0.0, 10.j, 10.0j, 0.0, 2.917-21.217j, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 10.0j, 10.0j, 0.0, 1.363-20.369j]], dtype = complex) 

# Step 2 Add model of load admittances

# Step 3 calculate internal gen voltages
# 2nd order model calc (following book example)
# Ei<delta = Vterm_mag + jXdi*(Pg-jQg)/Vterm 

#change Ei to vector , def. xd
#xd = 0.18
#Ei = np.add((Vmag[5],np.multiply(Qbus[2],np.divide(xd[2]/Vmag[5])))) #+ (Pbus[2]*xd[2]/Vmag[5])*1j   

def calcEi(Vmag,Vtheta,xd,P,Q):
    a = xd * Q 
    b = a / Vmag 
    Ei_r = Vmag + b #real part
    
    c = P * xd 
    Ei_i = c / Vmag #imag part
    
    Ei = Ei_r+Ei_i*1j
    #cmath.polar()Ei_r2+Ei_i2*1j
    Ei_mag = abs(Ei)
    phase_it = cmath.phase(Ei)  #angle between E and Vterm in rad
    gen_angle = phase_it + Vtheta #internal gen angle in rad
    return Ei_mag,gen_angle

Vmag[0], Vtheta[0] = calcEi(Vmag[4-1],Vtheta[4-1],Xd_trans[0],Pbus[1-1],Qbus[1-1]) 
Vmag[1], Vtheta[1] = calcEi(Vmag[5-1],Vtheta[5-1],Xd_trans[1],Pbus[2-1],Qbus[2-1])  
Vmag[2], Vtheta[2] = calcEi(Vmag[6-1],Vtheta[6-1],Xd_trans[2],Pbus[3-1],Qbus[3-1])

# Step 4 cacluate prefault / fault / post fault admittance matrices 

# Step 5 Kron reduction 
# Yreduced = Ynn - Yns*(1/Yss)*Ysn
# n = NGENS , s = remaining system nodes
def kronRed(Y,n,s): 
    Ynn = Y[0:n,0:n]
    Yns = Y[0:n,n:s] 
    Ysn = Y[n:s,0:n]
    Yss_inv = inv(Y[n:s,n:s])
    Yhat = Ynn - np.dot(np.dot(Yns,Yss_inv),Ysn)
    return Yhat 

Ypre_red = kronRed(Ybus_pre,NGEN,NBUS)
Yfault_red = kronRed(Ybus_fault,NGEN,NBUS)
Ypost_red = kronRed(Ybus_post,NGEN,NBUS)
#===============================================================================
# Solve Equations with scscipy.integrate.solve_ivp(f) --> NOT WORKING
# Solve Equations with Forward Euler
# M*delta(w') = Pm - Pe - D*delta(w)
# d' = delta(w)
# Ttransdo*Etransq' = Ef - Etransq + Id*(Xd - Xtransd)

# delta(w) = 'speed deviation' --> variable dOmega
#===============================================================================
#
#Solving process:
#For(whole time)
    # init w/ pre-fault condits
    # from 0 to 0.1 solve with fault matrix
    # from 0.1 to  1.5s solve with post fault matrix  

# #time points
timepoints = np.linspace(0,1.5)  #change time steps 
 
#Calc speed deviation
gen1 = np.zeros([np.size(timepoints),2])  #gen1 is an array of [rotor angles, speed deviatons] 

#Derivative Functions
def dOmega_hat(t,dWdt,Pm,Pe,D,dOmega,H):
    M = 0.6 #2*H*S/w_s 
    dWdt = 1/M*(Pm - Pe - D*dOmega) 
    return dWdt

def delta_hat(t,delta): #function prob not necessary for delta_hat = dOmega 
    dOmega = Omega[t]
    return dOmega 

def emfTransQ_hat(Ef,Eq_trans,Id,Xd,Xd_trans,Tdo_trans):
    return dEq_trans_dt

def emfTransD_hat(Ef,Ed_trans,Iq,Xq,Xq_trans,Tqo_trans):
    return dEq_trans_dt

def gen_Model(t,y,c):
    H = 10 #2*H*S/w_s 
    M = 0.6
    Pm = 1.999 
    Pe = 1.75
    D = 1.0 ## does this need to be here
    dWdt = 1/M*(Pm - Pe - D*y[0]) 
    derivs = [dWdt,y[0]]
    return derivs
# # # init condits
gen1[0,:] = Vtheta[0], 0.0  #init condits [delta0, w0] --> not used here
# gen2[0,:] = Vtheta[1], 0.0 
# gen3[0,:] = Vtheta[2], 0.0
c = 1
sol = solve_ivp(lambda t, y: gen_Model(t,y,c),[0,1.5],
                gen1[0,:],t_eval=timepoints)

#Solve for rotor angle
# sol = solve_ivp(delta_hat,[0, 1.5],Vtheta[0],t_eval=timepoints)
# sol_gen2 = solve_ivp(delta_hat,[0,1.5],Vtheta[1],t_eval=timepoints)
# sol_gen3 = solve_ivp(delta_hat,[0,1.5],Vtheta[2],t_eval=timepoints)

#  #plot rotor angles
plt.plot(sol.t,sol.y[0,:])
plt.plot(sol.t,sol.y[1,:])
# plt.plot(sol_gen2.t,sol_gen2.y[0,:])
# plt.plot(sol_gen3.t,sol_gen3.y[0,:])
#plt.plot(t,response[:,1])
plt.xlabel('time')
plt.ylabel('Rotor Angle')
plt.show()



# 
# def Pg_i (gen,Ybus,Vmag,Vtheta):
#     #np.real(Ybus[ii]) = Gii, np.imag(Ybus[ii]) = B
#     Pg = (Vmag[gen-1] ** 2) * np.real(Ybus[gen-1,gen-1])
#     # Pe = (Vmag[gen]^2 * G[ii]) + Vmag[gen]*Vmag[gen_k]*(B[ik]*sin(Vtheta[gen]-Vtheta[k])... sum
#     return Pg

#use this call to send arguments 
gen1_w = solve_ivp(lambda t, y: (t,dWdt,Pm,Pe,D,dOmega,H),[0,1.5],[0.0],t_eval=timepoints)

print(gen1_w)
# #forward Euler Yn+1 = Yn + h*f(Tn,Yn)
# h = 0.01 #step size
# 
# # #time points
# t = np.linspace(0,1.5)  #change time steps 
# 
# #Calc speed deviation
# gen1 = np.zeros([np.size(t),2])  #gen1 is an array of [rotor angles, speed deviatons] 
# 
# # # init condits
# gen1[0,:] = Vtheta[0], 0.0  #init condits [delta0, w0]
# 
# #Iterate through time steps
# for tstep in range(len(t)-1):
#     gen1[tstep+1,0] = gen1[tstep,0] + h*gen1[tstep,1]  #iterate next step for rotor angle calc (this is delta_hat at time t)
#     
#     gen1[tstep+1,1] = gen1[tstep,0] + h*dOmega_hat(tstep,Pm,Pe,D,dOmega) 
    




# ODE solver method -- NOT WORKING
##Gen Model --> expand each equation to be a vector representing each gen; or loop function for each gen
# def genModel (gens, t):
#     swingEqn, deltaW = gens  # vector of variable outputs
#     M = 0.6 #2*H*S/w_s 
#     Pm = 1.999 
#     Pe = 1.75
#     DAMP = 1.0 ## does this need to be here?
# #     swingEqn = (1/M)*Pm - Pe - DAMP*deltaW # M*delta(w') = Pm - Pe - D*delta(w)
# #     speedEqn = deltaW
#     dgdt = [(1/M)*Pm - Pe - DAMP*deltaW, deltaW] #, emfEqn]
#     return dgdt # dgdt = synch gen model diff eqs
# 
# def Pg_i (gen,Ybus,Vmag,Vtheta):
#     #np.real(Ybus[ii]) = Gii, np.imag(Ybus[ii]) = B
#     Pg = (Vmag[gen-1] ** 2) * np.real(Ybus[gen-1,gen-1])
#     # Pe = (Vmag[gen]^2 * G[ii]) + Vmag[gen]*Vmag[gen_k]*(B[ik]*sin(Vtheta[gen]-Vtheta[k])... sum
#     return Pg
# 
# Pg1 = Pg_i(1,Ypre_red,Vmag,Vtheta)
# 
# # init condits
# gen0 = [Vtheta[0], 0.0]  #init condits [delta0, w0]
# 
# #time points
# t = np.linspace(0,1.5)  #change time steps 
# 
# #solve gen eqns
# response = odeint(genModel,gen0,t)
#  
#  #plot
# plt.plot(t,response[:,0])
# plt.plot(t,response[:,1])
# plt.xlabel('time')
# plt.ylabel('y(t)')
# plt.show()
