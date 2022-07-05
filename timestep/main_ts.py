import numpy as np
import math 

g = 9.82665
timestep = 0.01
Sd = 58
TOW = 10
W = TOW*g
AD = 1000
rho0 = 1.225
rho = 1.225*(1+((-0.0065*AD)/288.15))**(-(1/(-0.0065))*((9.086/287.0531)-0.0065))
mi = 0.05

# Dados Aerodinâmicos da ASA
Sw = 2
iw = 3
iw = iw*math.pi/180
Clmaxw = 2.1
Cl0w = 0.2
Clalphaw = 0.05
Clalphaw = Clalphaw*180/math.pi
Cmacw = -0.08
MACw = 0.45

# Dados do Estabilizador horizontal
Sh = 0.3


# Distâncias (Em relação ao CG)
x_tdp = 0.40
z_tdp = 0.080
x_tdn = 0.160
x_t = 0.180
z_t = 0.250
x_ac = 0.020 # positivo para centro aerodinâmico na frente do cg


Rroda = 0.040

# Condições iniciais
i=0
t[i] = 0
x[i] = 0
z[i] = 0
vx[i] = 0
vz[i] = 0
V[i] = math.sqrt(vx[i]**2 + vz[i]**2)
teta[i] = 0
q[i] = 0
q_dot[i] =0


T[i] = 0
Lw[i] = 0.5*rho*V[i]**2*Sw*Clw
Lh[i] = 0.5*rho*V[i]**2*Sh*(Clh+de_cl_h) # tem que ver certinho 
R_n[i] =  -T[i]*np.sin(teta[i]) - L[i] + W
R_tdn[i] = ((x_tdp/(x_tdp+x_tdn))*W - Lw[i]*(x_tdn-x_ac) +Lh*()  )/x_tdn



while True:
    if R_n>=0: # Aeronave ainda em contato com o solo
        if R_tdn>0: # Aeronave em corrida
            pass
        else: # Aeronave rotacionando
    else: # Aeronave em subida





        if 



