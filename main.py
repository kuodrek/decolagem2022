import numpy as np
import math

# Descobrir o que é gama, teta e alfa nesse contexto

g = 9.81
rho = 1.086
m = 1
x_tdp = 1
x_tdn = 1
W = m*g
Iyy = 0.3
dt = 0
Sref = 1

def tracao(V, rho):
    pass

def takeoff(X, U, aircraft_data):
    Vx = X[0]
    Vz = X[1]
    q = X[2]
    x_pos = X[3]
    z_pos = X[4]
    teta = X[5]

    de = U[0]

    CL_alfa = aircraft_data['CL_alfa']
    CL_de = aircraft_data['CL_de']
    CL_q = aircraft_data['CL_q']
    CL_0 = aircraft_data['CL_0']
    Cm_alfa = aircraft_data['Cm_alfa']
    Cm_de = aircraft_data['Cm_de']
    Cm_q = aircraft_data['Cm_q']
    Cm_0 = aircraft_data['Cm_0']

    CL = CL_alfa * alfa + CL_de * de + CL_q * q + CL_0
    CD = 1
    Cm_cg = Cm_alfa * alfa + Cm_de * de + Cm_q * q + Cm_0

    V = math.sqrt(Vx ** 2 + Vz ** 2)

    T = tracao(V, rho) * g
    M_motor = T * dt

    L = 0.5 * rho * V ** 2 * Sref * CL
    D = 0.5 * rho * V ** 2 * Sref * CD
    M_cg = 0.5 * rho * V ** 2 * Sref * Cm_cg


    # Equações das reações dos trens de pouso / nariz
    # R_n = R_tdp + R_tdn
    R_n = W - T*np.sin(alfa) - L
    R_tdn = (R_n - 1 / x_tdp * (M_cg + M_motor)) * (1 + x_tdn/x_tdp) ** -1
    R_tdp = 1 /(x_tdp) * (M_cg + M_motor + R_tdn*x_tdn)
    P_cg = -R_tdp*x_tdp + R_tdn*x_tdn

    if R_n >= 0:
        # caso onde o avião não saiu do chão
        # Referencia: runway
        ax = (T*np.cos(alfa) - D + W - mi*R_n)/m
        az = (T*np.sin(alfa) + L - W + R_n)/m
        ateta = (M_cg + M_motor + P_cg) / Iyy
    else:
        # caso onde o avião saiu do chão
        # Referencia: horizon
        gama = teta - alfa
        ax = (T*np.cos(teta) - D*np.cos(gama) - L*np.sin(gama))/m
        az = (T*np.sin(teta) + L*np.cos(gama) - W - D*np.sin(gama))/m
        ateta = (M_cg + M_motor)/Iyy