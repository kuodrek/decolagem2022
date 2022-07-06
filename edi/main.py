from turtle import update
import numpy as np
from functions import takeoff, update_X
from data import get_data, get_data_teste

dados_planilha = get_data_teste()
de_takeoff = dados_planilha["de_takeoff"]

X_0 = np.zeros(6)
U_0 = [0]
U = U_0
h = 0.25

ac_eh = dados_planilha["ac_eh"]

t = 0
t_vetor = [t]
t_max = 20
X_vetor = [X_0]
data_output_vetor = []
X = X_0
while True:
    x_pos = X[3]
    if x_pos > ac_eh:
        U = [de_takeoff]

    K1, data_output = takeoff(X, U, dados_planilha)
    K2, _ = takeoff(X, U, dados_planilha)
    K3, _ = takeoff(X, U, dados_planilha)
    K4, _ = takeoff(X, U, dados_planilha)
    X = update_X(X, K1, K2, K3, K4, h)

    X_vetor.append(X)
    data_output_vetor.append(data_output)

    t += h
    t_vetor.append(t)
    if t > t_max:
        break
