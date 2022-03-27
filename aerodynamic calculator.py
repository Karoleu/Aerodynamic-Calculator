# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 20:07:55 2022

@author: Karol
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt
import math



D_rotor = 13.5 #srednica wirnika
Rtip = D_rotor*0.5
pH = 1 # [atm]
TH = 250
P_nominal = 670*1000


HubTipRatio = 0.1
rho = pH*101325/(287*TH) # gestosc atmosfery
u_tip = 265 # predkosc obwodowa wierzchołka lopaty
omega = u_tip/Rtip
TOW = 5000*9.81 # ciezar 
A = math.pi*0.25*D_rotor*D_rotor
ChordRoot = 119*0.001
ChordTip = 90*0.001
Theta0 = 11.17*math.pi/180
Nblades = 4
Cl_slope = 5.7
alfa_0 = -5*math.pi/180
SoundSpeed = math.sqrt(1.4*287*TH)#340
V_pozioma = 2/3.6
P_engines = P_nominal*pH/1.01325*math.sqrt(288.15/TH)
V_pionowa = 0
A_czolowe = 6
Cd_platowca = 0.4
Foporu = Cd_platowca*A_czolowe*0.5*rho*V_pozioma**2
i_angle = math.atan(Foporu/TOW) # kąt inklinacji wirnika w radianach
T = TOW
v_hover = math.sqrt(T/(2*rho*A))
v_AD = -0.5*V_pionowa + math.sqrt((0.5*V_pionowa)**2 + v_hover**2) # predkosc indukowana dla modelu actuator disk
P_AD = T*(V_pionowa+v_AD)
dynamic_viscosity = 1.8e-5
#print("v_hover = %1.5e"%v_hover)
#print("Sila oporu = %1.4e N\t tan_i = %1.3e\ti_angle = %2.3f stopni \t Area = %3.1f m^2"%(Foporu, Foporu/TOW, i_angle*180/math.pi, A))
n_stations = 501
BLADE = np.zeros((n_stations,19))
#zmienna: r/R r  C  Theta  u_P  u_T  fi  AoA  Cl  Cd  p_dyn L'  D',  T'  Q'  Re Mach  Pin'dMG
#indeks:   0  1  2    3    4    5    6   7    8    9  10    11  12   13  14  15  16  17    18
tabela = np.zeros((20,3))

Vmin = 10
Vmax = 100

EXCESS_THRUST = np.zeros((100,2)) 

def main():
    global V_pionowa
    
    Y = FindTheta(T)
    print(Y)
    print(" solved theta = %1.3f" % Y[0])
    Performance = Analyze(Y[0])
    print("performance: ")
    print(Performance)
    print("Pshaft/P_AD %1.4f"% (Performance[3]/P_AD))
    print( "# -----------------------------#")
    
    """
    Find_Vinduced(Chord, 4, Theta, 0.15)
    Ytip = Find_VinducedTip(Chord, Theta)
    print(Ytip)
    """   
   
    
    plt.plot(BLADE[:,0], BLADE[:,7]*180/math.pi)
    plt.plot(BLADE[:,0], BLADE[:,13])
    plt.yscale('log')
    plt.plot(EXCESS_THRUST[:,0],EXCESS_THRUST[:,1])
    plt.grid()
    
    


def FindTheta(RequiredThrust): # znajduje wartosc kąta theta dla zadanego ciagu
    def error(theta):
        VprimTip = Find_VinducedTip(ChordTip, theta)
        uPtip = V_pionowa + VprimTip[0]
        uTtip = omega*Rtip*0.9999 - VprimTip[2]
        sinFiT = math.sqrt((uPtip*uPtip)/(uPtip*uPtip+uTtip*uTtip))
        
        for i in range(len(BLADE)):
            BLADE[i,0] = HubTipRatio + (1-HubTipRatio)*i/(n_stations-1)
            BLADE[i,1] = BLADE[i,0] * Rtip
            BLADE[i,2] = ChordRoot + (ChordTip-ChordRoot)*i/(n_stations-1)
            BLADE[i,3] = theta
            vPrime = Find_Vinduced(BLADE[i,2], BLADE[i,1],BLADE[i,3], sinFiT)
            BLADE[i,4] = V_pionowa + vPrime[0]
            BLADE[i,5] = omega * BLADE[i,1] - vPrime[1]
            BLADE[i,6] = math.atan(BLADE[i,4]/BLADE[i,5])
            BLADE[i,7] = BLADE[i,3] - BLADE[i,6]
            BLADE[i,15] = BLADE[i,2]*math.sqrt(BLADE[i,4]*BLADE[i,4]+BLADE[i,5]*BLADE[i,5])*rho/dynamic_viscosity
            BLADE[i,16] = math.sqrt(BLADE[i,4]**2 + BLADE[i,5]**2)/SoundSpeed
            BLADE[i,8] = CL_od_alfa(BLADE[i,7], BLADE[i,16])
            BLADE[i,9] = CD_od_alfa(BLADE[i,7], BLADE[i,16])
            BLADE[i,10] = rho*0.5*(BLADE[i,4]*BLADE[i,4]+BLADE[i,5]*BLADE[i,5])
            BLADE[i,11] = BLADE[i,8]*BLADE[i,10]*BLADE[i,2]
            BLADE[i,12] = BLADE[i,9]*BLADE[i,10]*BLADE[i,2]
            if i == len(BLADE)-1:
                BLADE[i,13] = 0.0
            else:
                BLADE[i,13] = Nblades*(BLADE[i,11]*math.cos(BLADE[i,6]) - BLADE[i,12]*math.sin(BLADE[i,6])) 
            BLADE[i,14] = Nblades*BLADE[i,1]*(BLADE[i,11]*math.sin(BLADE[i,6]) + BLADE[i,12]*math.cos(BLADE[i,6]))
            BLADE[i,17] = BLADE[i,13] * BLADE[i,4]#Pinduced_prime
            BLADE[i,18] = BLADE[i,11]
            
        Thrust = np.trapz(BLADE[:,13], BLADE[:, 1])
        #Torque = np.trapz(BLADE[:,14], BLADE[:, 1])
        #Mgnacy = np.trapz(BLADE[:,18], BLADE[:, 1])
        #Pinduced = np.trapz(BLADE[:,17], BLADE[:, 1])
        #P_shaft = Torque * omega
        #print("theta = %1.3f, T = %1.3f"%(theta*180/math.pi, Thrust))
        return RequiredThrust - Thrust
    
    RESULT = opt.toms748(error, 0, 0.5, full_output = True)
    #print(RESULT)
    return RESULT[0], RESULT[1].converged
    
    

    
    """
    for j in range(10):
        MM = j/9*1.1
        DD = np.zeros((100, 2))
        for i in range(100):
            alfa = -6 + i/99*12
            alfa1 = alfa*math.pi/180
            DD[i,0] = alfa
            DD[i,1] = CD_od_alfa(alfa1, MM)
  
        plt.plot(DD[:,0], DD[:,1])
        #plt.yscale('log')
    plt.grid()
    
    DD = np.zeros((100, 2))
    for i in range(100):
        rR = 0.1 + 0.8999*i/99
        sinffiT = math.sin(10*math.pi/180)
        DD[i,0] = rR
        DD[i,1] = F_Prandtl(rR, Nblades, sinffiT)
    plt.plot(DD[:,0], DD[:,1])
    #plt.yscale('log')
    plt.grid()
    """
    
def Analyze(theta):

    VprimTip = Find_VinducedTip(ChordTip, theta)
    uPtip = V_pionowa + VprimTip[0]
    uTtip = omega*Rtip*0.9999 - VprimTip[2]
    sinFiT = math.sqrt((uPtip*uPtip)/(uPtip*uPtip+uTtip*uTtip))
    
    for i in range(len(BLADE)):
        BLADE[i,0] = HubTipRatio + (1-HubTipRatio)*i/(n_stations-1)
        BLADE[i,1] = BLADE[i,0] * Rtip
        BLADE[i,2] = ChordRoot + (ChordTip-ChordRoot)*i/(n_stations-1)
        BLADE[i,3] = theta
        vPrime = Find_Vinduced(BLADE[i,2], BLADE[i,1],BLADE[i,3], sinFiT)
        BLADE[i,4] = V_pionowa + vPrime[0]
        BLADE[i,5] = omega * BLADE[i,1] - vPrime[1]
        BLADE[i,6] = math.atan(BLADE[i,4]/BLADE[i,5])
        BLADE[i,7] = BLADE[i,3] - BLADE[i,6]
        BLADE[i,15] = BLADE[i,2]*math.sqrt(BLADE[i,4]*BLADE[i,4]+BLADE[i,5]*BLADE[i,5])*rho/dynamic_viscosity
        BLADE[i,16] = math.sqrt(BLADE[i,4]**2 + BLADE[i,5]**2)/SoundSpeed
        BLADE[i,8] = CL_od_alfa(BLADE[i,7], BLADE[i,16])
        BLADE[i,9] = CD_od_alfa(BLADE[i,7], BLADE[i,16])
        BLADE[i,10] = rho*0.5*(BLADE[i,4]*BLADE[i,4]+BLADE[i,5]*BLADE[i,5])
        BLADE[i,11] = BLADE[i,8]*BLADE[i,10]*BLADE[i,2]
        BLADE[i,12] = BLADE[i,9]*BLADE[i,10]*BLADE[i,2]
        if i == len(BLADE)-1:
            BLADE[i,13] = 0.0
        else:
            BLADE[i,13] = Nblades*(BLADE[i,11]*math.cos(BLADE[i,6]) - BLADE[i,12]*math.sin(BLADE[i,6])) 
        BLADE[i,14] = Nblades*BLADE[i,1]*(BLADE[i,11]*math.sin(BLADE[i,6]) + BLADE[i,12]*math.cos(BLADE[i,6]))
        BLADE[i,17] = BLADE[i,13] * BLADE[i,4]#Pinduced_prime
        BLADE[i,18] = BLADE[i,11]
        
    Thrust = np.trapz(BLADE[:,13], BLADE[:, 1])
    Torque = np.trapz(BLADE[:,14], BLADE[:, 1])
    #Mgnacy = np.trapz(BLADE[:,18], BLADE[:, 1])
    Pinduced = np.trapz(BLADE[:,17], BLADE[:, 1])
    P_shaft = Torque * omega
    print("TOW: %2.3f [kN], Thrust: %2.3f [kN], Pshaft = %1.1f [kW], Pinduced - %1.1f [kW], P_AD - %1.1f [kW], extraPower = %1.3f"  % (TOW/1000, Thrust/1000, P_shaft/1000, Pinduced/1000, P_AD/1000, (P_engines-P_shaft)/1000))
    return np.array([Thrust, Torque, Pinduced, P_shaft])

    
def Find_uP(Chord, r, Theta):
    
    def error(uP):
        fi = math.atan(uP/(omega*r))
        AoA  = Theta - fi
        CL = CL_od_alfa(AoA, 200000)
        CD = CD_od_alfa(AoA, 200000)
        q = rho*0.5 * (uP*uP + (omega*r)**2)
        Tprime = Nblades*q*Chord*(CL*math.cos(AoA)-CD*math.sin(AoA))
        if Tprime > 0 :
            uPresult = math.sqrt(Tprime/(4*rho*math.pi*r))
        elif Tprime < 0 :
            uPresult =  - math.sqrt(-Tprime/(4*rho*math.pi*r))
        else: 
            uPresult = 0.0
        return uP - uPresult
    v_h = math.sqrt(TOW/(2*rho*A))
    uP = opt.ridder(error, 0.0001, 6*v_h, full_output=True)
    #print("r: %1.3e"%r)
    #print(uP)
    
    return uP[0]

def Find_Vinduced(Chord, r, Theta, sinFiT):
    rR = r/Rtip #względny promień
    solidity = Nblades*Chord/(2*math.pi*r)
    FP = F_Prandtl(rR, Nblades, sinFiT)
    x0 = np.zeros(2)
    x0[0] = solidity*(V_pionowa**2 + (omega*r)**2)*0.6/(4*FP*v_hover)
    x0[1] = solidity*(V_pionowa**2 + (omega*r)**2)*0.015/(4*FP*v_hover)
    def error(X):
        vPprim = X[0]
        vTprim = X[1]
        uP = V_pionowa + vPprim
        uT = omega*r - vTprim
        U = math.sqrt(uT*uT+uP*uP)
        Mach = U/SoundSpeed
        sinFi = uP/U
        cosFi = uT/U
        fi = math.asin(uP/U)
        alfa = Theta - fi
        CL = CL_od_alfa(alfa, Mach)
        CD = CD_od_alfa(alfa, Mach)
        Cz = CL*cosFi - CD*sinFi
        Cx = CL*sinFi + CD*cosFi
        ER1 = solidity*U*U*Cz - 4*uP*vPprim*FP
        ER2 = solidity*U*U*Cx - 4*uP*vTprim*FP
        return np.array([ER1, ER2]) 
    
    RESULT = opt.fsolve(error, x0, full_output = True)
    if RESULT[2] != 1:
        print("Vinduced error\t r/R = %1.3f\t FP = %1.3f\tfsolve mes: %s"%(rR, FP, RESULT[3]))
    Y = [RESULT[0][0], RESULT[0][1], RESULT[2]]
    return Y

def Find_VinducedTip(Chord, Theta):
    rR = 0.9999
    r = rR*Rtip
    solidity = Nblades*Chord/(2*math.pi*r)
    x0 = np.zeros(2)
    x0[0] = solidity*(V_pionowa**2 + (omega*r)**2)*0.6/(4*0.1*v_hover)
    x0[1] = solidity*(V_pionowa**2 + (omega*r)**2)*0.015/(4*0.1*v_hover)
    def error(X):
        vPprim = X[0]
        vTprim = X[1]
        uP = V_pionowa + vPprim
        uT = omega*r - vTprim
        U = math.sqrt(uT*uT+uP*uP)
        Mach = U/SoundSpeed
        sinFi = uP/U
        cosFi = uT/U
        fi = math.asin(uP/U)
        alfa = Theta - fi
        CL = CL_od_alfa(alfa, Mach)
        CD = CD_od_alfa(alfa, Mach)
        FP = F_Prandtl(rR, Nblades, sinFi)
        Cz = CL*cosFi - CD*sinFi
        Cx = CL*sinFi + CD*cosFi
        ER1 = solidity*U*U*Cz - 4*uP*vPprim*FP
        ER2 = solidity*U*U*Cx - 4*uP*vTprim*FP
        return np.array([ER1, ER2]) 
    
    RESULT = opt.fsolve(error, x0, full_output = True)
    if RESULT[2] != 1:
        print("Vinduced error\t r/R = %1.3f\t fsolve mes: %s"%(rR, RESULT[3]))
    Y = [RESULT[0][0], RESULT[0][1], RESULT[2]]
    return Y
        
        

def F_Prandtl(r_R, N, sin_fi_T):
    if sin_fi_T < 0.01:
        sin_fi_T = 0.01
    if r_R >= 1.0:
        #print(" GIven r/R is higher than 1.0, r/R = %1.3e"%r_R)
        r_R =0.9999
    if r_R <= 0:
        #print(" GIven r/R is higher than 1.0, r/R = %1.3e"%r_R)
        r_R =0.0001
    ff = 0.5*N*(1-r_R)/sin_fi_T
    FP = 2/math.pi*math.acos(math.exp(-ff))
    return FP
 
def CL_od_alfa(AoA, MM):
    if MM < 0.2:
        M = 0.2
    elif MM > 0.75:
        M = 0.75
    else:
        M = MM
    
    def W4(B, A1, A2, A4, X):
        return  B+ A1*X + A2*X*X + A4*X*X*X*X
    
    alfa = AoA*180/math.pi
    alfa_crit = W4(19.04767426,	-31.27605946,	17.10937007,	-10.65943521, M)
    CL_crit = W4(1.864199647,	-2.726075386,	1.575964056,	-0.720560707, M)
    alfa_0 = W4(-1.069771767,	0.246802366,	0.060646154,	-0.01436723, M)
    
    Clinear = CL_crit*(alfa-alfa_0)/(alfa_crit-alfa_0)
    if alfa <= alfa_crit:
        CL = Clinear
    else:
        A = W4(-3.084463857,	4.501025852,	9.655165084,	-17.49634242, M)
        B = W4(-0.279478875,	-0.14237066,	0.069390413,	-0.073347099, M)
        C = W4( 3.664935018,	-0.146953544,	-1.230426913,	-0.612097531 ,M)
        D = W4( 5.556945648,	-0.23504624,	-0.325071786,	-0.625171243 ,M)
        x = alfa - alfa_crit
        CL = Clinear + (A*x**2 + B*x**3)/(D + C*x**2)
    return CL
    
def CD_od_alfa(AoA, M):
    CD_0 = 0.008
    alfa = abs(AoA*180/math.pi)
    deltaCD = 0
    if alfa > 10 :
        alfa = 10
    if alfa > 4 :
        x = alfa - 4
        deltaCD = 0.001*x + 0.002*x*x
    return CD_0 + deltaCD








"""
for i in range(101):
    V_pozioma = Vmin + i/100*(Vmax-Vmin)
    Foporu = Cd_platowca*A_czolowe*0.5*rho*V_pozioma**2
    i_angle = math.atan(Foporu/TOW)
    T = math.sqrt(Foporu**2 + TOW**2)
    v_h = math.sqrt(T/(2*rho*A))
    X = Find_v(v_h, V_pozioma, i_angle)
    print(X)
    tabela[i,0] = V_pozioma
    tabela[i,1:3] = X
    #tabela[i,2] = X[1]
"""    

"""
plt.plot(tabela[:,0], tabela[:,1])
plt.grid()
"""
    
main()