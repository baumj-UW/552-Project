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
Xd = np.array([[0.08, 0.016],
             [0.18, 0.036 ],
             [0.12, 0.024]])  #Array of synchronous reactances, [Xd, X'd] per gen
Xq = np.array([[0.16, 0.16],
             [0.36, 0.36 ],
             [0.24, 0.24]])  #Array of synchronous reactances, [Xq, X'q] per gen
#Tdo = np.array([[1.2],[1.4],[1.5]]) #Open circuit time constants T'do
Tdo = np.array([[5.0],[5.0],[5.0]]) #Open circuit time transient d time constants T'do
Tqo = np.array([[10.0],[10.0],[10.0]]) #Open circuit subtransient q time constants T"qo
D = np.array([[0.0],[0.0],[0.0]]) #Open circuit time constants T'do

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
               Ybus_pre_B[bus_i,bus_j] = -1/Xd[bus_i,0]  
               Ybus_pre_B[bus_j,bus_i] = -1/Xd[bus_i,0]
           else:
               Ybus_pre_B[bus_i,bus_j] = 1/Xd[bus_i,0]
               Ybus_pre_B[bus_j,bus_i] = 1/Xd[bus_i,0]

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
 
#"remove" bus 7 during fault <--- change this calc to use X'd!!!
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

def calcEi(Vt,Vtheta,xd,P,Q): #returns Eab based on input condits
    Ei_ab = (Vt + xd*Q/Vt) + 1j*(P*xd / Vt)
    return Ei_ab #send complex calc of E back 

# Calculate internal gen voltages <-- this can be made more generic based on BUS_CONN
Eab = np.zeros((NGEN,1),dtype=complex) 

Eab[0] = calcEi(Vmag[4-1],Vtheta[4-1,0],Xd[0,0],Pbus[1-1],Qbus[1-1]) 
Eab[1] = calcEi(Vmag[5-1],Vtheta[5-1,0],Xd[1,0],Pbus[2-1],Qbus[2-1])  
Eab[2] = calcEi(Vmag[6-1],Vtheta[6-1,0],Xd[2,0],Pbus[3-1],Qbus[3-1])

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


def abTdq(delta): #returns (a,b) (d,q) transform matrix
    T = np.array([[-math.sin(delta),math.cos(delta)],
                 [math.cos(delta),math.sin(delta)]])
    return T


# removing Eab calc w/ Xf --> use for 4th order model?
# #Initial guess of Xf and X'q (x'd given)
#Xq_trans = 0.5*Xd_trans
#Xf = 0.5*(Xd_trans + Xq_trans) 

#Estimate values for Ed_f,Eq_f --> this should be done with Xf from above
# Start w/ Ed_f = Eq_f = Ef from fault # guessing Edq = Eab w/o rotor angle adjustment..?
Ef = np.zeros((NGEN,1)) #Ef for E'd calc, constant set to init Eab magnitude 
Edq_f = np.zeros((NGEN,2)) #Array of [Ed_f, Eq_f], rows correspond to gens
Eab_v = np.zeros((NGEN,2)) #array of [Ea, Eb], rows correspond to gens
deltaT = np.zeros((NGEN,1)) #array of rotor angles to use in transform 
for gen in range(NGEN):
    Ef[gen] = abs(Eab[gen])
    deltaT[gen] = cmath.phase(Eab[gen]) + Vtheta[gen+NGEN,0] #actual init rotor angle 
    Eab_v[gen,:] = np.array([[np.real(Eab[gen,0]), np.imag(Eab[gen,0])]]) #Eab = Eab_v 
    Edq_f[gen,:] = np.dot(abTdq(deltaT[gen]),Eab_v[gen,:].T).T

#Calculate generator currents
 #constant load admittance, expect Edq as a row [E'd, E'q] per gen 
 # Edq is the TRANSIENT emf 
 # deltaT is array of rotor angles for dq transform
def Ig(Ybus,Edq,deltaT):   #takes E'dq and returns Idq 
    Ig_ab = np.zeros((NGEN,2))
    Ig_dq = np.zeros((NGEN,2))
    for gen_i in range(NGEN):
        # dq to ab for current calc on common axis
        Ea,Eb = np.dot(abTdq(deltaT[gen_i]),Edq[gen_i,:].T)
        for gen_j in range(NGEN):
            Gij = np.real(Ybus[gen_i,gen_j])
            Bij = np.imag(Ybus[gen_i,gen_j])
            Ig_ab[gen_i,0] += Gij*Ea - Bij*Eb
            Ig_ab[gen_i,1] += Gij*Eb + Bij*Ea 
        
        #convert Iab to Idq
        Ig_dq[gen_i,:] = np.dot(abTdq(deltaT[gen_i]),Ig_ab[gen_i,:].T)
        # 
        
    return Ig_dq #returns (NGENx2) array of [Id,Iq] per gen


# #Iteration to solve Edq process <-- skipping this for 3rd order
# Ig_dq_l = Ig(Yfault_red,Edq_f,deltaT) #Idq at iteration l
#  
# #Correct the Edq values for next iteration step 
# #Edq(l+1) = E'dq - (Xq_trans - Xf)
# for gen in range(NGEN):
#     mismatch_d = (Xq_trans[gen]-Xf[gen])*Ig_dq_l[gen,1]
#     mismatch_q = (Xd_trans[gen]-Xf[gen])*Ig_dq_l[gen,0]
#     Edq_f[gen,0] += mismatch_d
#     Edq_f[gen,1] -= mismatch_q
#      

# #time points
fault_times = np.linspace(0,F_CLEAR,round(F_CLEAR/0.005))  
postf_times = np.linspace(F_CLEAR,END_SIM,round((END_SIM-F_CLEAR)/0.005))

#remove Vmag from inputs, not needed for 3rd order
def gen_Model(t,y,Vmag,Ybus,Pm,M,Ef,Xd,Xq,Tdo,Tqo): #y is an array of state variables [w,delta,E'q,E'd]
    
    D = 0.0 ## Neglect Pd for 2nd order Model

    omega = y[0:NGEN]
    delta = y[NGEN:2*NGEN]
    EqTran = y[2*NGEN:3*NGEN]
    EdTran = y[3*NGEN:4*NGEN]
    
    #Use Edq to find Idq
    Edq = np.array((np.array(EdTran),np.array(EqTran))).T
    Emag = np.linalg.norm(Edq,axis=1)
    Idq = Ig(Ybus,Edq,delta.T) #returns [Id, Iq] for all gens
    
    Pe = np.zeros((NGEN,1))
    dWdt = np.zeros((NGEN,1))
    dEqTran = np.zeros((NGEN,1))
    dEdTran = np.zeros((NGEN,1))
    for gen in range(NGEN):
        Pe[gen] = Pg_i(gen,Ybus,Emag,delta)
        dWdt[gen] = 1/M[gen]*(Pm[gen] - Pe[gen] - D*omega[gen])  
        dEqTran[gen] = 1/Tdo[gen]*(Ef[gen] - EqTran[gen] + Idq[gen,0]*(Xd[gen,0] - Xd[gen,1])) # Xd columns [Xd, X'd]
        dEdTran[gen] = 1/Tqo[gen]*(-EdTran[gen] - Idq[gen,1]*(Xq[gen,0] - Xq[gen,1])) # Xd columns [Xd, X'd]
    derivs = [dWdt[0],dWdt[1],dWdt[2],omega[0],omega[1],omega[2],
              dEqTran[0],dEqTran[1],dEqTran[2],
              dEdTran[0],dEdTran[1],dEdTran[2]] 
    return derivs


def Pg_i (gen_i,Ybus,Vmag,delta): #assume gens are numbered {0,1,2}
    Pg = 0.0    #init Pg to 0 then calc.
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
initGen = np.zeros((2,4*NGEN)) #[w1, w2, w3, delta1, delta2, delta3, EqTran, EdTran..] w = 0
for gen in range(NGEN):
    initGen[0,gen+NGEN]  = deltaT[gen,0] 
    initGen[0,gen+2*NGEN] = Edq_f[gen,1] # 'init values of EqTran...' from internal emf calcs?
    initGen[0,gen+3*NGEN] = Edq_f[gen,0] # init values of EdTran for 4th order
    
#solve response during fault time 0,0.1s
fault_sol = solve_ivp(lambda t, y: gen_Model(t,y,abs(Eab),Yfault_red,Pbus[0:NGEN],M,Ef,Xd,Xq,Tdo,Tqo),
                [0,F_CLEAR],initGen[0,:],t_eval=fault_times)   

## Set new init condits from fault solution 
for gen in range(NGEN):
    initGen[1,gen] = fault_sol.y[gen,len(fault_sol.t)-1] #omega at F_CLEAR
    initGen[1,gen+NGEN] = fault_sol.y[gen+NGEN,len(fault_sol.t)-1] #delta at F_CLEAR
    initGen[1,gen+2*NGEN] = fault_sol.y[gen+2*NGEN,len(fault_sol.t)-1] #EqTran at F_CLEAR
    initGen[1,gen+3*NGEN] = fault_sol.y[gen+3*NGEN,len(fault_sol.t)-1] #EdTran at F_CLEAR

#solve response post fault clearing
postf_sol = solve_ivp(lambda t, y: gen_Model(t,y,abs(Eab),Ypost_red,Pbus[0:NGEN],M,Ef,Xd,Xq,Tdo,Tqo),
                [F_CLEAR,END_SIM],initGen[1,:],t_eval=postf_times)    




# Combine solution results to plot <-- remove repeated time step at fault clear
sim_times = np.concatenate((fault_sol.t,postf_sol.t))
results = np.zeros((4*NGEN,len(sim_times))) #array of results [speed;delta;Eq]

speedFig = plt.figure(1) #rotor speed plot
for omega in range(NGEN):
    results[omega,:] = np.concatenate((fault_sol.y[omega,:],
                                       postf_sol.y[omega,:]),axis=0)
    plt.plot(sim_times,(results[omega,:]+W_S)/W_S,
             label='Gen '+str(omega+1))    
plt.xlabel('Time (sec)')
plt.ylabel('Rotor Speed (pu)')
plt.legend()
plt.grid(True)

rotorAngleFig = plt.figure(2) #rotor angle plot
for delta in range(NGEN,NGEN*2):
    results[delta,:] = np.concatenate((fault_sol.y[delta,:],postf_sol.y[delta,:]))
    plt.plot(sim_times,(180/math.pi)*results[delta,:],
             label='Gen '+str(delta-NGEN+1))
plt.xlabel('Time (sec)')
plt.ylabel('Rotor Angle (deg)')
plt.legend()
plt.grid(True)


#Create relative rotor angle plot
relativeAngleFig = plt.figure(3)
rel_delta = np.zeros((NGEN-1,len(sim_times)))
for delta in range(NGEN-1):
    rel_delta[delta,:] = results[delta+NGEN+1,:] - results[NGEN,:] #relative angle (Gen_i - Gen1)
    plt.plot(sim_times,(180/math.pi)*rel_delta[delta,:],label='$\delta$'+str(delta+2)+'1')
plt.xlabel('Time (sec)')
plt.ylabel('Relative Rotor Angles (deg)')
plt.legend()
plt.grid(True)

emfFig = plt.figure(4)
#Create Eq plot
for Eq in range(NGEN*2,NGEN*3):
    results[Eq,:] = np.concatenate((fault_sol.y[Eq,:],postf_sol.y[Eq,:]))
    plt.subplot(2,1,1) #4th subplot in figs
    plt.plot(sim_times,results[Eq,:],label='Gen '+str(Eq-NGEN*2+1)) # add label to plots
plt.xlabel('Time (sec)')
plt.ylabel("E'q (pu)")
plt.legend()
plt.grid(True)

#Create Eq plot
for Ed in range(NGEN*3,NGEN*4):
    results[Ed,:] = np.concatenate((fault_sol.y[Ed,:],postf_sol.y[Ed,:]))
    plt.subplot(2,1,2) #4th subplot in figs
    plt.plot(sim_times,results[Ed,:],label='Gen '+str(Ed-NGEN*3+1)) # add label to plots
plt.xlabel('Time (sec)')
plt.ylabel("E'd (pu)")
plt.legend()
plt.grid(True)

plt.show()
    





