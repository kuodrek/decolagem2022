import numpy as np
from functions import takeoff
from data import get_data

dados_planilha = get_data()
de_takeoff = dados_planilha['de_takeoff']

X_0 = np.zeros(6)
U_0 = [0]
U = U_0
h = 0.25

x_acionamento = dados_planilha['x_acionamento']

t = 0
t_vetor = [t]
t_max = 20
while True:
    x_pos = X[3]
    if x_pos > x_acionamento:
        U = [de_takeoff]

    K1, data_output = takeoff(X, U, dados_planilha)
    K2, _ = takeoff(X, U, dados_planilha)
    K3, _ = takeoff(X, U, dados_planilha)
    K4, _ = takeoff(X, U, dados_planilha)
    X = X+h/6*(K1+2*K2+2*K3+K4)

    t += h
    t_vetor.append(t)
    if t > t_max: break
