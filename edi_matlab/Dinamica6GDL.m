function [Xp,CF,reacoes,Fa] = Dinamica6GDL(X,U,H,AircraftData,estado_do_aviao)
%% Dados do avi?o
geral = AircraftData{2,1};
controle = AircraftData{1,1};
g = geral(1,2);
m = controle(1,5);
W = m*g;
Ib = AircraftData{7,1};
Iyy = Ib(2,2);
x_tdp = geral(1,13);
x_tdn = geral(1,14);
mi = geral(1,15);
%% Vetor de estados (S? ? utilizado u, w, q, teta, x, z)
u = X(1);
v = X(2);
w = X(3);
% p = X(4);
q = X(5);
% r = X(6);
phi = X(10);
teta = X(11);
psi = X(12);
%% Velocidades e ?ngulos de ataque / derrapagem
alfa_max = 20 * pi / 180;
gama_max = 20 * pi / 180;
if strcmp(estado_do_aviao, 'subida')
    gama = atan(w/u);
    %% Restri??o de gama m?ximo
    if gama > gama_max
        gama = gama_max;
    end
    if gama < -gama_max
        gama = -gama_max;
    end
    %% Restri??o de alfa m?ximo
    alfa = teta - gama;
    if alfa > alfa_max
        alfa = alfa_max;
    end
    if alfa < -alfa_max
        alfa = -alfa_max;
    end
else
    alfa = teta;
end
%% For?as e momentos aerodin?micos e de tra??o
[Fa,Ft,Ha,Ht,CF] = fmsolver_decolagem2022(X,U,AircraftData,estado_do_aviao);

% For?as aerodin?micas
D = Fa(1);
L = Fa(3);

% For?a tra??o
T = Ft(1);

% Momentos aerodin?mico e de tra??o
Ma = Ha(2);
Mt = Ht(2);
%% Rea??o dos trens de pouso e nariz
R_n = W - T*sin(teta) - L;

if R_n < 0 || strcmp(estado_do_aviao, 'subida') == true
    R_n = 0;
end

R_tdn = (R_n - 1 / x_tdp * (Ma - Mt)) / (1 + x_tdn/x_tdp);
if R_tdn < 0 || strcmp(estado_do_aviao, 'rotacao') == true || strcmp(estado_do_aviao, 'subida') == true
    R_tdn = 0;
end

% R_tdp = 1 /(x_tdp) * (Ma - Mt + R_tdn*x_tdn);
R_tdp = R_n - R_tdn;
if R_tdp < 0 || strcmp(estado_do_aviao, 'subida') == true
    R_tdp = 0;
end
P_cg = -R_tdp*x_tdp + R_tdn*x_tdn;

if R_n > 0
    % Rea??o dos trens de pouso e nariz > 0 -> Avi?o correndo na pista
    % Referencia: runway
    up = (T* cos(teta) - D - mi*R_n)/m;
    wp = 0;
    if R_tdn > 0
        % Rea??o do trem de nariz > 0 -> Avi?o correndo na pista (n?o rotacionou)
        qp = 0;
        tetap = 0;
        x_pos_p = u;
        z_pos_p = 0;
    else
        % Rea??o do trem de nariz = 0 -> Avi?o rotacionando
        qp = (Ma - Mt + P_cg) / Iyy;
        tetap = q;
        x_pos_p = u;
        z_pos_p = 0;
    end
else
    % Rea??es dos trens de pouso e nariz = 0 -> Avi?o voando
    % Referencia: horizon
    
    gama = atan(w/u);
    alfa = teta - gama;
    if gama > gama_max
        gama = gama_max;
    end
    if gama < -gama_max
        gama = -gama_max;
    end
    if alfa > alfa_max
        alfa = alfa_max;
    end
    if alfa < -alfa_max
        alfa = -alfa_max;
    end
    
    up = (T* cos(teta) - D* cos(gama) - L* sin(gama))/m;
    wp = (T* sin(teta) + L* cos(gama) - W - D* sin(gama))/m;
    qp = (Ma - Mt)/Iyy;
    x_pos_p = u;
    z_pos_p = w;
    tetap = q;
end
%% Vetor de sa?da
reacoes = [R_n R_tdp R_tdn];
Xp = [up 0 wp 0 qp 0 x_pos_p 0 z_pos_p 0 tetap 0];
end