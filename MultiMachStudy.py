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
import numpy as np  # Methods for linear algebra
from numpy.linalg import inv
from scipy.integrate import odeint  #refs odeint directly instead of long pointer
from scipy.integrate import solve_ivp #ODE45 equiv (use this instead of odeint)
import matplotlib.pyplot as plt  #refs this pointer as plt --> try simplifiying this later


#Define Constants

NGEN = 3 
NBUS = 8
NLOAD = 2 #not currently used
F_CLEAR = 0.1 #Time at which the fault is cleared 
END_SIM = 1.5 #End time of simulation 
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

W_S = 2*math.pi*60 #synchronous speed 

Sn = 1.0 #Per unit rating of the generators (100MVA base) 
H = np.array([[10.0],[3.01],[6.4]]) # Inertia constant from book (100MVA base)
Xd_trans = np.array([[0.08],[0.18],[0.12]]) #transient reactance from book table

# Step 1 - Convert to common base (already done for this project) 

#Calculate Inertia Coefficient M
M = np.zeros((NGEN,1))
for gen in range(NGEN):
    M[gen] = 2*H[gen]*Sn/W_S

# Initialize empty bus informaiton arrays <-- split Vmag/theta into V and E?
Vmag = np.zeros((NBUS,1))   #vector of voltage magnitudes
Vtheta = np.zeros((NBUS,2))  # array of volt angles in degrees [theta @ t=0, theta @ t=F_CLEAR] <-- not needed
Pbus = np.zeros((NBUS,1)) #vector of P at each bus
Qbus = np.zeros((NBUS,1)) #vector of Q at each bus


#assign outputs from book solution  <-- replace w/ NR solution 
Vmag[NGEN:NBUS] = np.array([[1.04],[1.02],[1.05],[0.9911],[1.0135]])
## Degrees Vtheta[NGEN:NBUS] = [[0.0],[-3.55],[-2.90],[-7.48],[-7.05]]
Vtheta[NGEN:NBUS,0:1] = np.array([[0.0],[-0.06196],[-0.05061],[-0.13055],[-0.12305]])

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
            Ybus_pre_B[bus_i,bus_j] = 1/LINE_XL

#System bus self admittances             
for bus_i in range(NGEN,NBUS):
    Ybus_pre_B[bus_i,bus_i] = -1*np.sum(Ybus_pre_B[bus_i,:]) + LINE_YB*(np.sum(BUS_CONN[bus_i,NGEN:NBUS])-1)


            
#Calculate equiv gen bus E and delta
#Create full augmented Ybus matrix (includes gen and load admittance)
# Pre-fault: calc. from steady state augmented matrix
# fault: Assume 3ph-g fault, set 

### Step 2 Add model of load admittances 
#Calculate constant load bus admittances and append Ybus
# Load P and Q from solution <-- this will come from NR 
Pbus[6] = 2.8653
Qbus[6] = 1.2244

Pbus[7] = 1.4
Qbus[7] = 0.4

for load in range(6,8): # more generic way to iterate through loads?
    Ybus_pre_G[load,load] += Pbus[load]/(Vmag[load] ** 2) 
    Ybus_pre_B[load,load] -= Qbus[load]/(Vmag[load] ** 2)
    
#Create complex pre-fault Ybus    
Ybus_pre = Ybus_pre_G + 1j*Ybus_pre_B
 
#"remove" bus 7 during fault 
Ybus_fault = np.delete(Ybus_pre,6,0) #remove row 7 from Ybus
Ybus_fault = np.delete(Ybus_fault,6,1) #remove column 7 from Ybus

# Alter Ybus matrix post-fault, remove line 6-7
Ybus_post = Ybus_pre
Ybus_post[5,5] -= 1j*(LINE_YB - 1/LINE_XL) 
Ybus_post[6,6] -= 1j*(LINE_YB - 1/LINE_XL) #
Ybus_post[5,6] = 0.0
Ybus_post[6,5] = 0.0


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

# Calculate internal gen voltages <-- this can be made more generic based on BUS_CONN
Vmag[0], Vtheta[0,0] = calcEi(Vmag[4-1],Vtheta[4-1,0],Xd_trans[0],Pbus[1-1],Qbus[1-1]) 
Vmag[1], Vtheta[1,0] = calcEi(Vmag[5-1],Vtheta[5-1,0],Xd_trans[1],Pbus[2-1],Qbus[2-1])  
Vmag[2], Vtheta[2,0] = calcEi(Vmag[6-1],Vtheta[6-1,0],Xd_trans[2],Pbus[3-1],Qbus[3-1])

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
fault_times = np.linspace(0,F_CLEAR)  #change time steps 
postf_times = np.linspace(F_CLEAR,END_SIM)

##Derivative Functions
# def dOmega_hat(t,dWdt,Pm,Pe,D,dOmega,H):
#     M = 0.6 #2*H*S/w_s 
#     dWdt = 1/M*(Pm - Pe - D*dOmega) 
#     return dWdt
# 
# def delta_hat(t,delta): #function prob not necessary for delta_hat = dOmega 
#     dOmega = Omega[t]
#     return dOmega 
# 
# def emfTransQ_hat(Ef,Eq_trans,Id,Xd,Xd_trans,Tdo_trans):
#     return dEq_trans_dt
# 
# def emfTransD_hat(Ef,Ed_trans,Iq,Xq,Xq_trans,Tqo_trans):
#     return dEq_trans_dt

def gen_Model(t,y,Vmag,Ybus,Pm,M): #y is an array of [w1,w2,w3,d1,d2,d3]
    #H = 10 #2*H*S/w_s 
    #M = 0.6
    D = 0 ## Neglect Pd for 2nd order Model
    #Pm = 1.999 #<-- Pm from Pbus init state 

    omega = y[0:NGEN]
    delta = y[NGEN:]
    
    Pe = np.zeros((NGEN,1))
    dWdt = np.zeros((NGEN,1))
    for gen in range(NGEN):
        Pe[gen] = Pg_i(gen,Ybus,Vmag,delta)  
        dWdt[gen] = 1/M[gen]*(Pm[gen] - Pe[gen] - D*omega[gen])  #
 
    derivs = [dWdt[0],dWdt[1],dWdt[2],omega[0],omega[1],omega[2]] #simplify with concat[dWdt; omegas]
    return derivs





def Pg_i (gen_i,Ybus,Vmag,delta): #assume gens are numbered {0,1,2}, only send Vmag[gens] (create separate Emag array)
    Pg = 0.0    #init Pg to 0 then calc.
    #np.real(Ybus[ii]) = Gii, np.imag(Ybus[ii]) = B
    for gen_j in range(NGEN):
        Gii = np.real(Ybus[gen_i,gen_i])
        Gij = np.real(Ybus[gen_i,gen_j])
        Bij = np.imag(Ybus[gen_i,gen_j])
        di = delta[gen_i]
        dj = delta[gen_j]
        if gen_j == (gen_i):
            Pg += (Vmag[gen_i] ** 2) * Gii
        else:
            Pg += Vmag[gen_i]*Vmag[gen_j]*(Bij*math.sin(di-dj) + Gij*math.cos(di-dj))

    return Pg


# # Pre-fault initial condits and args --> send as an array to gen_Model 
#Create array of initial generator conditions [pre-fault;postfault clear]
initGen = np.zeros((2,2*NGEN)) #[w1, w2, w3, delta1, delta2, delta3] w = 0
for gen in range(NGEN):
    initGen[0,gen+NGEN]  = Vtheta[gen,0] 

#solve response during fault time 0,0.1s
fault_sol = solve_ivp(lambda t, y: gen_Model(t,y,Vmag[0:NGEN],Yfault_red,Pbus[0:NGEN],M),
                [0,F_CLEAR],initGen[0,:],t_eval=fault_times)   

## Set new init condits from fault solution 
for gen in range(NGEN):
    initGen[1,gen] = fault_sol.y[gen,len(fault_sol.t)-1] #omega at F_CLEAR
    initGen[1,gen+NGEN] = fault_sol.y[gen+NGEN,len(fault_sol.t)-1] #delta at F_CLEAR

#solve response post fault clearing
postf_sol = solve_ivp(lambda t, y: gen_Model(t,y,Vmag[0:NGEN],Ypost_red,Pbus[0:NGEN],M),
                [F_CLEAR,END_SIM],initGen[1,:],t_eval=postf_times)    



figs = plt.figure(1)
# Combine solution results to plot <-- remove repeated time step at fault clear
sim_times = np.concatenate((fault_sol.t,postf_sol.t))
results = np.zeros((2*NGEN,len(sim_times))) #array of results [speed;delta]
for omega in range(NGEN):
    results[omega,:] = np.concatenate((fault_sol.y[omega,:],
                                       postf_sol.y[omega,:]),axis=0)
    plt.subplot(2,2,1) #first subplot in figs
    plt.plot(sim_times,(results[omega,:]+W_S)/W_S,
             label='Gen '+str(omega+1)) # add label to plots    
plt.xlabel('Time (sec)')
plt.ylabel('Rotor Speed (pu)')
plt.legend()
plt.grid(True)

for delta in range(NGEN,NGEN*2):
    results[delta,:] = np.concatenate((fault_sol.y[delta,:],postf_sol.y[delta,:]))
    plt.subplot(2,2,2) #2nd subplot in figs
    plt.plot(sim_times,(180/math.pi)*results[delta,:],
             label='Gen '+str(delta-NGEN+1)) # add label to plots
plt.xlabel('Time (sec)')
plt.ylabel('Rotor Angle (deg)')
plt.legend()
plt.grid(True)

#Create relative rotor angle plot
rel_delta = np.zeros((NGEN-1,len(sim_times)))
for delta in range(NGEN-1):
    rel_delta[delta,:] = results[delta+NGEN+1,:] - results[NGEN,:] #relative angle (Gen_i - Gen1)
    plt.subplot(2,2,3) #3rd subplot in figs
    plt.plot(sim_times,(180/math.pi)*rel_delta[delta,:],label='$\delta$'+str(delta+2)+'1')
plt.xlabel('Time (sec)')
plt.ylabel('Relative Rotor Angles (deg)')
plt.legend()
plt.grid(True)

plt.show(figs)
    

# cmath.rect(r, phi) --> to convert polar to rect. value to combine Vmag and Vtheta

## plot concats
# #  #plot speeds
# plt.plot(sol.t,sol.y[0,:],label='Gen1 w')
# plt.plot(sol.t,sol.y[1,:],label='Gen2 w')
# plt.plot(sol.t,sol.y[2,:],label='Gen3 w')
# # plt.plot(sol_gen2.t,sol_gen2.y[0,:])
# # plt.plot(sol_gen3.t,sol_gen3.y[0,:])
# #plt.plot(t,response[:,1])
# plt.xlabel('time')
# plt.ylabel('Speed')
# plt.legend()
# plt.show()

#  #plot speeds
#plt.plot(sol.t,sol.y[0,:],label='Gen1 w')
# 


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
