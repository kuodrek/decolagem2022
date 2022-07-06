import numpy as np
import math

# Descobrir o que é gama, teta e alfa nesse contexto

def update_X(X,K1,K2,K3,K4,h):
    X_temp = np.zeros([len(X)])
    for i, _ in enumerate(X):
        X_temp[i]=X[i]+h/6 * (K1[i]+2*K2[i]+2*K3[i]+K4[i])
    return X_temp

def tracao(V, rho, helice_dados):
    rho0 = 1.225
    T = rho / rho0 * np.polyval(helice_dados,V)
    return T

def takeoff(X, U, dados_planilha):
    # X = [Vx, Vz, q, x_pos, z_pos, teta]
    # X_dot = [ax, az, ateta, Vx, Vz, q]
    # U = [de]

    Vx = X[0]
    Vz = X[1]
    q = X[2]
    x_pos = X[3]
    z_pos = X[4]
    teta = X[5]

    de = U[0]

    if Vz == 0: 
        alfa = teta
        gama = 0
    else: 
        alfa = np.arctan(Vz/Vx)
        gama = teta - alfa

    V = math.sqrt(Vx ** 2 + Vz ** 2)

    # Dados avião
    g=dados_planilha['g']
    rho=dados_planilha['rho']
    m=dados_planilha['m']
    x_tdp=dados_planilha['x_tdp']
    x_tdn=dados_planilha['x_tdn']
    Iyy=dados_planilha['Iyy']
    dt=dados_planilha['dt']
    Sref=dados_planilha['Sref']
    mi=dados_planilha['mi']
    W=m*g

    CL_alfa = dados_planilha['CL_alfa']
    CL_de = dados_planilha['CL_de']
    CL_q = dados_planilha['CL_q']
    CL_0 = dados_planilha['CL_0']
    Cm_alfa = dados_planilha['Cm_alfa']
    Cm_de = dados_planilha['Cm_de']
    Cm_q = dados_planilha['Cm_q']
    Cm_0 = dados_planilha['Cm_0']
    helice_dados = dados_planilha['helice_dados']
    CD_poly = dados_planilha['CD_poly']

    # Coeficientes
    CL = CL_alfa * alfa + CL_de * de + CL_q * q + CL_0
    CD = np.polyval(CD_poly, alfa)
    Cm_cg = Cm_alfa * alfa + Cm_de * de + Cm_q * q + Cm_0

    # Força e momento do motor
    T = tracao(V, rho, helice_dados) * g
    M_motor = T * dt

    # Forças e momentos aerodinâmicos
    L = 0.5 * rho * V ** 2 * Sref * CL
    D = 0.5 * rho * V ** 2 * Sref * CD
    M_cg = 0.5 * rho * V ** 2 * Sref * Cm_cg

    # Equações das reações dos trens de pouso / nariz
    # R_n = R_tdp + R_tdn
    R_n = W - T*np.sin(teta) - L
    R_tdn = (R_n - 1 / x_tdp * (M_cg + M_motor)) * (1 + x_tdn/x_tdp) ** -1
    if R_tdn < 0: R_tdn = 0

    R_tdp = 1 /(x_tdp) * (M_cg + M_motor + R_tdn*x_tdn)
    P_cg = -R_tdp*x_tdp + R_tdn*x_tdn

    if R_n >= 0:
        # Reação dos trens de pouso e nariz > 0 -> Avião correndo na pista
        # Referencia: runway
        ax = (T*np.cos(teta) - D - mi*R_n)/m
        az = 0
        if R_tdn > 0:
            # Reação do trem de nariz > 0 -> Avião correndo na pista (não rotacionou)
            ateta = 0
            teta_dot = 0
        else:
            # Reação do trem de nariz = 0 -> Avião rotacionando
            ateta = (M_cg + M_motor + P_cg) / Iyy
            teta_dot = q
        x_pos_dot = Vx
        z_pos_dot = 0
    else:
        # Reações dos trens de pouso e nariz = 0 -> Avião voando
        # Referencia: horizon
        ax = (T*np.cos(teta) - D*np.cos(gama) - L*np.sin(gama))/m
        az = (T*np.sin(teta) + L*np.cos(gama) - W - D*np.sin(gama))/m
        ateta = (M_cg + M_motor)/Iyy
        x_pos_dot = Vx
        z_pos_dot = Vz
        teta_dot = q
    
    X_dot = np.zeros(6)
    X_dot[0] = ax
    X_dot[1] = az
    X_dot[2] = ateta
    X_dot[3] = x_pos_dot
    X_dot[4] = z_pos_dot
    X_dot[5] = teta_dot

    dados_output = {
        'R_n': R_n,
        'R_tdp': R_tdp,
        'R_tdn': R_tdn,
        'gama': gama,
        'alfa': alfa,
        'L': L,
        'D': D,
        'T': T,
        'M_cg': M_cg
    }

    return X_dot, dados_output
