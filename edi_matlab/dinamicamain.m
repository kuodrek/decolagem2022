% DINAMICA 6 GDL SETUP
clear
clc
clf
%% Chamar datasetup para os dados da aeronave
disp('Chamando datasetup')
global data_check
global AircraftData
global SetupVoo
if isempty(data_check) == 1
    data_check = 1;
    [data,setup] = datasetup();
    AircraftData = data;
    SetupVoo = setup;
end
disp('Dados adquiridos!')
%% Vetores iniciais X_0 e U_0
X_0 = zeros(1,12);
superficies = AircraftData{3,1};
dtmax = superficies(4,6);
U_0 = [dtmax 0 0 0];
H_0 = [0 0 0];
% Vetor de reações dos trens de pouso / nariz iniciais. Precisam de um
% valor arbitrario pra iniciar a simulação. Logo após a 2a iteração é
% calculado o valor correto
reacoes = [999 999 999];
% A variável estado_do_aviao pode assumir três valores:
% 'corrida' -> 1a parte da decolagem
% 'rotacao' -> 2a parte da decolagem
% 'subida' -> 3a parte da decolagem
estado_do_aviao = 'corrida';
%% Solver numérico
disp('INICIANDO SIMULAÇÃO')
t_inicial = 0;
t_final = 10;
h=20e-4;
n_pto = (t_final-t_inicial)/h+1;
vet_t=linspace(t_inicial,t_final,n_pto);   % Vetor dos tempos amostrados
U_vet=zeros(n_pto,4);
vetor_de_estados={};
solucao = zeros(n_pto,12);
Xp_vetor = zeros(n_pto,12);
X = X_0;
U = U_0;
H = H_0;
check1 = 0;
check2 = 0;

Sd = 58;
ac_eh = 43; 
de_takeoff = 0 * pi / 180;

ac_desativacao = 53;

for i=1:n_pto
    solucao(i,:) = X;
    t = vet_t(i);
    
    %% Verificação da fase da decolagem
    R_n = reacoes(1);
    R_tdp = reacoes(2);
    R_tdn = reacoes(3);
    
    x_pos = X(7);
    if x_pos > ac_eh
        U(2) = de_takeoff;
    end
    
    if x_pos > ac_desativacao
        U(2) = 0;
    end
    
    if R_tdn > 0 && R_tdp > 0
        estado_do_aviao = 'corrida';
    elseif R_tdn == 0 && R_tdp > 0
        estado_do_aviao = 'rotacao';
    elseif R_tdn == 0 && R_tdp == 0
        estado_do_aviao = 'subida';
    end
    vetor_de_estados{i} = {estado_do_aviao};
%% Runge-kutta
    U_vet(i,:) = U;
    K1=Dinamica6GDL(X,U,H,AircraftData,estado_do_aviao);
    [~,CF,reacoes,Fa]=Dinamica6GDL(X,U,H,AircraftData,estado_do_aviao);
    CF_vetor(i,:) = CF;
    K2=Dinamica6GDL(X+h/2*K1,U,H,AircraftData,estado_do_aviao);
    K3=Dinamica6GDL(X+h/2*K2,U,H,AircraftData,estado_do_aviao);
    K4=Dinamica6GDL(X+h*K3,U,H,AircraftData,estado_do_aviao);
    X=X+h/6*(K1+2*K2+2*K3+K4);
    
    L_vetor(i,:) = Fa(3);
    reacoes_vetor(i,:) = [reacoes x_pos];
    Xp_vetor(i,:) = K1;
    %% Condição de parada
%     R_n = reacoes(1);
%     if R_n < 0
%         break
%     end
end
%% Pós processamento
t_excel = vet_t';
disp('Fim da simulação. Plotando gráficos...')
%  Xp = [up vp wp pp qp rp xp_e yp_e zp_e phip tetap psip]';
u = solucao(:,1);
v = solucao(:,2);
w = solucao(:,3);
V = sqrt(u.^2+v.^2+w.^2);
alfa = atan(w./u)*180/pi;
beta = atan(v./V)*180/pi;
p = solucao(:,4)*180/pi;
q = solucao(:,5)*180/pi;
r = solucao(:,6)*180/pi;
x = solucao(:,7);
y = solucao(:,8);
z = solucao(:,9);
phi = solucao(:,10)*180/pi;
teta = solucao(:,11)*180/pi;
psi = solucao(:,12)*180/pi;

figure(1)
subplot(2,3,1);
plot(vet_t,reacoes_vetor(:,1));
title('Reacao total x t')
xlim([0 t_final])
ylim([0 max(reacoes_vetor(:,1))])

subplot(2,3,2);
plot(vet_t,reacoes_vetor(:,2));
title('Reacao tdp x t')
xlim([0 t_final])
ylim([0 max(reacoes_vetor(:,2))])

subplot(2,3,3);
plot(vet_t,reacoes_vetor(:,3));
title('Reacao tdn x t')
xlim([0 t_final])
ylim([0 max(reacoes_vetor(:,3))])

subplot(2,3,4);
plot(x,z);
title('Z x X')
xlim([0 max(x)])
ylim([0 40])

subplot(2,3,5);
teta_0 = teta(1)*ones(size(teta,1),size(teta,2));
plot(vet_t,teta);
title('teta x t')
xlim([0 t_final])
ylim([-30 30])

subplot(2,3,6);
q_0 = q(1)*ones(size(q,1),size(q,2));
plot(vet_t,q);
title('q x t')
xlim([0 t_final])
ylim([-20 20])

figure(2)
plot(vet_t,L_vetor);
title('L x t')
xlim([0 t_final])
ylim([0 1.2*max(L_vetor)])
