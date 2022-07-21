% DECOLAGEM MONOPLANO 2022 REGULAR
clear % Deleta vari�veis do workspace
clc % Limpa command window
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

massa_inicial = 13;
massa_final = 15;
massa_step = 0.25;
massa_vetor = (massa_inicial:massa_step:massa_final)';
% massa_vetor = [14.4];

de_inicial = -15;
de_final = -10;
de_step = 1;
de_vetor = (de_inicial:de_step:de_final)';

acionamento_inicial = 30;
acionamento_final = 38;
acionamento_step = 2;
acionamento_vetor = (acionamento_inicial:acionamento_step:acionamento_final)';
output = cell(size(acionamento_vetor,1),1);
z_threshold = 0.8;
%% Pegar valores de itera��o
for ii=1:size(de_vetor,1)
    de_takeoff =  de_vetor(ii) * pi / 180;
    for jj=1:size(acionamento_vetor,1)
        ac_eh = acionamento_vetor(jj);
        massa_max = 0;
        for kk=1:size(massa_vetor,1)
            %% Varia��o da massa
            m = massa_vetor(kk);
            AircraftData{1,1}(1,5) = m;
            %% Vari�veis de deflex�o do profundor
            x_inicial = ac_eh;
            x_final = ac_eh;
            %% Vetores iniciais X_0 e U_0
            X_0 = zeros(1,12);
            geral = AircraftData{2,1};
            z_inicial = geral(1,17) ; % Altura do cg inicial[m]
            X_0(9) = z_inicial;
            superficies = AircraftData{3,1};
            dtmax = superficies(4,6);
            U_0 = [dtmax 0 0 0];
            H_0 = [0 0 0];
            % Vetor de rea��es dos trens de pouso / nariz iniciais. Precisam de um
            % valor arbitrario pra iniciar a simula��o. Logo ap�s a 2a itera��o �
            % calculado o valor correto
            reacoes = [999 999 999];
            % A vari�vel estado_do_aviao pode assumir tr�s valores:
            % 'corrida' -> 1a parte da decolagem
            % 'rotacao' -> 2a parte da decolagem
            % 'subida' -> 3a parte da decolagem
            estado_do_aviao = 'corrida';
            %% Valores gerais
            superficies = AircraftData{3,1};
            iw = superficies(1,7);
            %% Solver num�rico
            t_inicial = 0;
            t_final = 20;
            h=20e-3;
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
            %% Defini��o de valores limite
            % Velocidade de arfagem limite
            qmax = 15* pi / 180;
            %% Vari�veis de decolagem ao liftoff (Avi�o sai do ch�o)
            x_liftoff = 0;
            V_liftoff = 0;
            gama_liftoff = 0;
            teta_liftoff = 0;
            q_liftof = 0;
            % Vari�vel que pega o estado do avi�o ao passar pelo obst�culo.
            % liftoff_check = false -> Avi�o n�o saiu do ch�o
            % liftoff_check = true -> Avi�o saiu do ch�o
            liftoff_check = false;
            %% Vari�veis de decolagem ao avi�o passar pelo obst�culo
            x_decolage=0;
            Vz_decolage=0;
            V_decolage=0;
            z_decolage=0;
            gama_decolage=0;
            teta_decolage=0;
            alfa_decolage=0;
            % Vari�vel que pega o estado do avi�o ao passar pelo obst�culo.
            % obstaculo_check = false -> Avi�o n�o passou pelo obst�culo
            % obstaculo_check = true -> Avi�o passou pelo obst�culo
            obstaculo_check = false;
            %% Simula��o
            for i=1:n_pto
                %% Gravar resultados da itera��o atual
                t = vet_t(i);
                solucao(i,:) = X;
                x_pos = X(7);
                %% Ativa��o/Desativa��o do profundor
                if x_pos > ac_eh % fuguinha == 1 -> profundor � ativado
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
                % Faz sentido esse crit�rio final?
                % Desativar o profundor ap�s o final da decolagem
                if x_pos > Sd
                    U(2) = -2*pi/180;
                end
                %% Verifica��o do estado do avi�o
                R_tdp = reacoes(2);
                R_tdn = reacoes(3);
                if R_tdn > 0 && R_tdp > 0
                    estado_do_aviao = 'corrida';
                elseif R_tdn == 0 && R_tdp > 0
                    estado_do_aviao = 'rotacao';
                elseif R_tdn == 0 && R_tdp == 0
                    estado_do_aviao = 'subida';
                    %% Pegar estado do avi�o no momento em que ele sai do ch�o
                    if liftoff_check == false
                        x_liftoff = x_pos;
                        V_liftoff = sqrt(X(1).^2+X(3).^2);
                        teta_liftoff = X(11);
                        q_liftoff = X(5);
                        liftoff_check = true;
                    end
                    %% Verificar se o avi�o passa ou n�o do obst�culo em 58m (Altura do obst�culo: 0.7m)
                    if abs(x_pos-Sd)<0.05 && obstaculo_check == false
                        x_decolage = x_pos;
                        Vz_decolage = X(3);
                        V_decolage = sqrt(X(1).^2+X(3).^2);
                        z_decolage = X(9);
                        gama_decolage = atan(X(3)/X(1));
                        teta_decolage = X(11);
                        alfa_decolage = teta_decolage - gama_decolage;
                        obstaculo_check = true;
                    end
                end
                %% Integrador Runge-kutta
                U_vet(i,:) = U;
                K1=Dinamica6GDL(X,U,H,AircraftData,estado_do_aviao);
                [~,~,reacoes,~]=Dinamica6GDL(X,U,H,AircraftData,estado_do_aviao);
                % Vetor de coeficientes aerodin�micos
                K2=Dinamica6GDL(X+h/2*K1,U,H,AircraftData,estado_do_aviao);
                K3=Dinamica6GDL(X+h/2*K2,U,H,AircraftData,estado_do_aviao);
                K4=Dinamica6GDL(X+h*K3,U,H,AircraftData,estado_do_aviao);
                X=X+h/6*(K1+2*K2+2*K3+K4);
                %% Restri��es f�sicas
                % Verificar se |q| > |qmax| no estado atual do avi�o
                if X(5) > qmax
                    X(5) = qmax;
                end
                if X(5) < -qmax
                    X(5) = -qmax;
                end
                
                % Se o avi�o est� correndo ele n�o pode ter teta < 0
                if X(11) < 0 && (strcmp(estado_do_aviao, 'corrida') && strcmp(estado_do_aviao, 'rotacao'))
                    X(11) = 0;
                end
                %% Vetores de for�as e momentos
                reacoes_vetor(i,:) = [reacoes x_pos];
            end
            fprintf("\n================================");
            fprintf("\nRodando para ac_eh: %gm",ac_eh);
            fprintf("\nRodando para massa: %gkg",m);
            fprintf("\nRodando para deflexao: %g�",de_takeoff * 180 / pi);
            fprintf("\nz_decolage: %gm", z_decolage);
            fprintf("\n================================");
            if z_decolage > z_threshold && m > massa_max
                output{ii,1}(jj,:) = [ac_eh m z_decolage];
                output{ii,2} = de_takeoff * 180 / pi;
                massa_max = m;
            end
        end
    end
end

%% Salvar dados no excel
salvar_dados(output)