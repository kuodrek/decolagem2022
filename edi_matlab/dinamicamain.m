% DECOLAGEM MONOPLANO 2022 REGULAR
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
t_final = 15;
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
Sd = geral(1,16);
%% Definição de valores limite
% Velocidade de arfagem limite
qmax = 10* pi / 180;
%% Variáveis de deflexão do profundor
ac_eh = 30;
de_takeoff = -10* pi / 180;
x_inicial = ac_eh;
x_final = ac_eh+5;
%% Variáveis de decolagem ao liftoff (Avião sai do chão)
x_liftoff = 0;
V_liftoff = 0;
gama_liftoff = 0;
teta_liftoff = 0;
q_liftof = 0;
% Variável que pega o estado do avião ao passar pelo obstáculo.
% liftoff_check = false -> Avião não saiu do chão
% liftoff_check = true -> Avião saiu do chão
liftoff_check = false;
%% Variáveis de decolagem ao avião passar pelo obstáculo
x_decolage=0;
Vz_decolage=0;
V_decolage=0;
z_decolage=0;
gama_decolage=0;
teta_decolage=0;
alfa_decolage=0;
% Variável que pega o estado do avião ao passar pelo obstáculo.
% obstaculo_check = false -> Avião não passou pelo obstáculo
% obstaculo_check = true -> Avião passou pelo obstáculo
obstaculo_check = false;
%% Simulação
for i=1:n_pto
    %% Gravar resultados da iteração atual
    t = vet_t(i);
    solucao(i,:) = X;
    x_pos = X(7);
    %% Ativação/Desativação do profundor
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
    if x_pos > Sd
        U(2) = -2*pi/180;
    end
    %% Verificação do estado do avião
    R_tdp = reacoes(2);
    R_tdn = reacoes(3);
    if R_tdn > 0 && R_tdp > 0
        estado_do_aviao = 'corrida';
    elseif R_tdn == 0 && R_tdp > 0
        estado_do_aviao = 'rotacao';
    elseif R_tdn == 0 && R_tdp == 0
        estado_do_aviao = 'subida';
        %% Pegar estado do avião no momento em que ele sai do chão
        if liftoff_check == false
            x_liftoff = x_pos;
            V_liftoff = sqrt(X(1).^2+X(3).^2);
            teta_liftoff = X(11);
            q_liftoff = X(5);
            liftoff_check = true;
        end
        %% Verificar se o avião passa ou não do obstáculo em 58m (Altura do obstáculo: 0.7m)
        if abs(x_pos-Sd)<0.05 && obstaculo_check == false
            x_decolage = x_pos;
            Vz_decolage = X(3);
            V_decolage = sqrt(X(1).^2+X(3).^2);
            z_decolage = X(9);
            gama_decolage = atan(X(3)/X(1));
            teta_decolage = X(11);
            alfa_decolage = teta_decolage - gama_decolage;
            if z_decolage < 0.7
                fprintf('gol do marcilio dias :(\n')
            end
            obstaculo_check = true;
        end
    end
    vetor_de_estados{i,1} = {estado_do_aviao};
    vetor_de_estados{i,2} = t;
    %% Integrador Runge-kutta
    U_vet(i,:) = U;
    K1=Dinamica6GDL(X,U,H,AircraftData,estado_do_aviao);
    [~,CF,reacoes,Fa]=Dinamica6GDL(X,U,H,AircraftData,estado_do_aviao);
    % Vetor de coeficientes aerodinâmicos
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
    %% Vetores de forças e momentos
    L_vetor(i,:) = Fa(3);
    reacoes_vetor(i,:) = [reacoes x_pos];
    Xp_vetor(i,:) = K1;
end
%% Pós processamento
t_excel = vet_t';
disp('Fim da simulação. Plotando gráficos...')
%  Xp = [up vp wp pp qp rp xp_e yp_e zp_e phip tetap psip]';
u = solucao(:,1);
w = solucao(:,3);
V = sqrt(u.^2+w.^2);
gama = atan(w./u)*180/pi;
q = solucao(:,5)*180/pi;
x = solucao(:,7);
z = solucao(:,9);
teta = solucao(:,11)*180/pi;
alfa = teta - gama;
%% Output da decolagem
fprintf('---LIFTOFF---\n')
alfa_liftoff = teta_liftoff - gama_liftoff;
fprintf('x[m] ao decolar: %.2f\n', x_liftoff)
fprintf('V[m/s] ao decolar: %.2f\n', V_liftoff)
fprintf('teta[°] ao decolar: %.2f\n', teta_liftoff * 180/pi)
fprintf('alfa da asa[°] ao decolar: %.2f\n', (teta_liftoff+iw) * 180/pi)
fprintf('q[°/s] ao decolar: %.2f\n', q_liftoff * 180 / pi)
fprintf('---OBSTÁCULO---\n')
fprintf('z[m] aos %gm: %.2f\n', Sd, z_decolage)
fprintf('x[m] aos %gm: %.2f\n', Sd, x_decolage)
fprintf('Vz[m/s] aos %gm: %.2f\n', Sd, Vz_decolage)
fprintf('V[m/s] aos %gm: %.2f\n', Sd, V_decolage)
fprintf('gama[°] aos %gm: %.2f\n', Sd, gama_decolage * 180 / pi)
fprintf('teta[°] aos %gm: %.2f\n', Sd, teta_decolage * 180 / pi)
fprintf('alfa da asa[°] aos %gm: %.2f\n', Sd, (alfa_decolage+iw) * 180 / pi)
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

subplot(2,1,2);
plot(x,z);
title('Z x X')
xlim([0 max(x)])
ylim([0 3*max(z)])

% subplot(2,3,3);
% teta_0 = teta(1)*ones(size(teta,1),size(teta,2));
% plot(x,teta);
% title('teta x X')
% xlim([0 max(x)])
% ylim([0 20])
%
% subplot(2,3,4);
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
hold on
plot(x,gama);
plot(x,U_vet(:,2)*180/pi)
hold off
title('gama x X')
legend('gama','\delta_{e}')
xlim([0 max(x)])
ylim([-30 30])

alfa_asa = alfa + ones(size(alfa,1),size(alfa,2))*iw*180/pi;
subplot(2,3,5)
hold on
plot(x,alfa);
plot(x,alfa_asa);
hold off
legend('alfa_{aeronave}','alfa_{asa}')
title('alfa x X')
xlim([0 max(x)])
ylim([-30 30])