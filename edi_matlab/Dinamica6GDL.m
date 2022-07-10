function [Xp,CF,reacoes,Fa] = Dinamica6GDL(X,U,H,AircraftData,estado_do_aviao)
%% Dados do avião
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
%% Vetor de estados (Só é utilizado u, w, q, teta, x, z)
u = X(1);
v = X(2);
w = X(3);
% p = X(4);
q = X(5);
% r = X(6);
phi = X(10);
teta = X(11);
psi = X(12);

%% Vento no sistema NED (north, east, down)
Nvento = H(1);
Evento = H(2);
Dvento = H(3);
% Matriz de transformação do sistema terrestre móvel para o sistema do corpo (aula 4 de ec nos slides do ITA)
Cbv = [cos(teta)*cos(psi),cos(teta)*sin(psi),-sin(teta);...
    sin(phi)*cos(teta)*cos(psi)-cos(phi)*sin(psi),cos(phi)*cos(psi)+sin(phi)*sin(teta)*sin(psi),sin(phi)*cos(teta);...
    cos(phi)*sin(teta)*cos(psi)+sin(phi)*sin(psi),-sin(phi)*cos(psi)+cos(phi)*sin(teta)*sin(psi),cos(phi)*cos(teta)];
% Vento no sistema do corpo da aeronave (X, Y, Z)
Bvento = Cbv*[Nvento;Evento;Dvento];
uvento = Bvento(1);
vvento = Bvento(2);
wvento = Bvento(3);

u = u + uvento;
v = v + vvento;
w = w + wvento;
X(1) = u;
X(2) = v;
X(3) = w;

%% Velocidades e ângulos de ataque / derrapagem
beta = 0;
alfa_max = 20 * pi / 180;
gama_max = 20 * pi / 180;
if strcmp(estado_do_aviao, 'subida')
    gama = atan(w/u);
    %% Restrição de gama máximo
    if gama > gama_max
        gama = gama_max;
    end
    if gama < -gama_max
        gama = -gama_max;
    end
    %% Restrição de alfa máximo
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
%% Forças e momentos aerodinâmicos e de tração
[Fa,Ft,Ha,Ht,CF] = fmsolver_decolagem2022(X,U,AircraftData,estado_do_aviao);

% Forças aerodinâmicas
D = Fa(1);
C = Fa(2);
L = Fa(3);

% Força tração
T = Ft(1);

% Momentos aerodinâmico e de tração
Ma = Ha(2);
Mt = Ht(2);
%% Reação dos trens de pouso e nariz
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
    % Reação dos trens de pouso e nariz > 0 -> Avião correndo na pista
    % Referencia: runway
    up = (T* cos(teta) - D - mi*R_n)/m;
    wp = 0;
    if R_tdn > 0
        % Reação do trem de nariz > 0 -> Avião correndo na pista (não rotacionou)
        qp = 0;
        tetap = 0;
        x_pos_p = u;
        z_pos_p = 0;
    else
        % Reação do trem de nariz = 0 -> Avião rotacionando
        qp = (Ma - Mt + P_cg) / Iyy;
        tetap = q;
        x_pos_p = u;
        z_pos_p = 0;
    end
else
    % Reações dos trens de pouso e nariz = 0 -> Avião voando
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

%% Acelerações da dinamica normal
% up = 1/m*(Xa+Xt) - w*q -g*sin(teta);
% wp = 1/m*(Za+Zt) + u*q +g*cos(teta);
% qp = 1/Iy*(Ma-Mt);
% np = u*cos(teta)+w*sin(teta);
% hp = u*sin(teta)-w*cos(teta);
% tetap = q;

%% Vetor de saída
if qp < 0
%     fprintf('Ma: %g\n', Ma)
%    fprintf('Mt: %g\n', Mt)
    aaaa = 1;
end
reacoes = [R_n R_tdp R_tdn];
Xp = [up 0 wp 0 qp 0 x_pos_p 0 z_pos_p 0 tetap 0];
end