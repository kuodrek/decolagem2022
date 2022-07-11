% DECOLAGEM MONOPLANO 2022 MICRO
clear % Deleta variáveis do workspace
clc % Limpa command window
clf('reset') % Reseta todos os gráficos
close all % Fecha todas as janelas de gráficos
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
geral = AircraftData{2,1};
z_inicial = geral(1,17) ; % Altura do cg inicial[m]
X_0(9) = z_inicial;
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
%% Valores gerais
superficies = AircraftData{3,1};
iw = superficies(1,7);
%% Solver numérico
disp('INICIANDO SIMULAÇÃO')
t_inicial = 0;
t_final = 5;
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
Sd = geral(1,16);

z_warning = 0.2;
%% Definição de valores de limite
qmax = 10 * pi / 180;
gama_max = 20 * pi / 180;
alfa_max = 20 * pi / 180;
%% Variáveis de deflexão do profundor
ac_eh = Sd;
de_takeoff = -14* pi / 180;
x_inicial = Sd;
x_final = Sd;
%% Simulação
for i=1:n_pto
    %% Gravar resultados da iteração atual
    t = vet_t(i);
    solucao(i,:) = X;
    x_pos = X(7);
    %% Verificação da fase da decolagem
    R_n = reacoes(1);
    R_tdp = reacoes(2);
    R_tdn = reacoes(3);
    if x_pos > ac_eh % fuguinha == 1 -> profundor é ativado
        if (x_final-x_inicial) == 0
            de_atual = de_takeoff;
        else
            de_atual = de_takeoff / (x_final - x_inicial) * (x_pos - x_inicial);
        end
        if x_pos > x_final
            de_atual = de_takeoff;
        end
        U(2) = de_atual;
    end
    % Final da decolagem: x_pos > Sd
    % Faz sentido esse critério final?
    % Desativar o profundor após o final da decolagem
    if x_pos > 10
        U(2) = -2*pi/180;
    end
    %% Verificação do estado do avião
    R_tdp = reacoes(2);
    R_tdn = reacoes(3);
    if R_tdn > 0 && R_tdp > 0
        estado_do_aviao = 'corrida';
    elseif R_tdn == 0 && R_tdp > 0
        estado_do_aviao = 'rotacao';
    elseif (R_tdn == 0 && R_tdp == 0) || x_pos > Sd
        estado_do_aviao = 'subida';
    end
    vetor_de_estados{i,1} = {estado_do_aviao};
    vetor_de_estados{i,2} = t;
%% Runge-kutta
    U_vet(i,:) = U;
    K1=Dinamica6GDL(X,U,H,AircraftData,estado_do_aviao);
    [~,CF,reacoes,Fa]=Dinamica6GDL(X,U,H,AircraftData,estado_do_aviao);
    CF_vetor(i,:) = CF;
    K2=Dinamica6GDL(X+h/2*K1,U,H,AircraftData,estado_do_aviao);
    K3=Dinamica6GDL(X+h/2*K2,U,H,AircraftData,estado_do_aviao);
    K4=Dinamica6GDL(X+h*K3,U,H,AircraftData,estado_do_aviao);
    X=X+h/6*(K1+2*K2+2*K3+K4);
    %% Restrições físicas
    % Verificar se |q| > |qmax| no estado atual do avião
    if X(5) > qmax
        X(5) = qmax;
    end
    if X(5) < -qmax
        X(5) = -qmax;
    end

    % Se o avião está correndo ele não pode ter teta < 0
    if X(11) < 0 && (strcmp(estado_do_aviao, 'corrida') && strcmp(estado_do_aviao, 'rotacao'))
        X(11) = 0;
    end
    %% Vetores de diversas forças e momentos
    L_vetor(i,:) = Fa(3);
    reacoes_vetor(i,:) = [reacoes x_pos t];
    Xp_vetor(i,:) = K1;
end
%% Pós processamento
t_excel = vet_t';
disp('Fim da simulação. Plotando gráficos...')
%  Xp = [up vp wp pp qp rp xp_e yp_e zp_e phip tetap psip]';
u = solucao(:,1);
v = solucao(:,2);
w = solucao(:,3);
V = sqrt(u.^2+v.^2+w.^2);
gama = atan(w./u)*180/pi;
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
alfa = teta - gama;
%% Gráficos
figure(1)
subplot(2,1,1);
hold on
plot(x,reacoes_vetor(:,1));
plot(x,reacoes_vetor(:,2));
plot(x,reacoes_vetor(:,3));
hold off
title('Reacões x X')
legend('R_{total}','R_{tdp}','R_{tdn}')
xlim([0 max(x)])
ylim([0 max(reacoes_vetor(:,1))])

z_warning_vetor = ones(size(z,1),size(z,2))*z_warning;
fprintf("zmin: %g\n", min(z))
subplot(2,1,2);
hold on
plot(x,z);
plot(x,z_warning_vetor);
title('Z x X')
xlim([0 30])
ylim([0 3])
hold off

% subplot(2,3,5);
% hold on
% teta_0 = teta(1)*ones(size(teta,1),size(teta,2));
% plot(x,teta);
% plot(x,U_vet(:,2)*180/pi)
% title('teta x X')
% xlim([0 30])
% ylim([-30 30])
% hold off
% 
% subplot(2,3,6);
% q_0 = q(1)*ones(size(q,1),size(q,2));
% plot(vet_t,q);
% title('q x t')
% xlim([0 t_final])
% ylim([-20 20])

figure(2)
subplot(2,3,1)
plot(vet_t,L_vetor);
title('L x t')
xlim([0 t_final])
ylim([0 1.2*max(L_vetor)])

subplot(2,3,2)
plot(vet_t,CF_vetor(:,3));
title('CL x t')
xlim([0 t_final])
ylim([0 3])

subplot(2,3,3)
plot(vet_t,V(:,1));
title('V x t')
xlim([0 t_final])
ylim([0 1.2*max(V)])

subplot(2,3,4)
plot(vet_t,gama);
title('gama x t')
xlim([0 t_final])
ylim([-30 30])

alfa_asa = alfa + ones(size(alfa,1),size(alfa,2))*iw*180/pi;
subplot(2,3,5)
hold on
plot(vet_t,alfa);
plot(vet_t,alfa_asa);
title('alfa x t')
xlim([0 t_final])
ylim([-30 30])
hold off