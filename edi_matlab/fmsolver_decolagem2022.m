function [Fa,Ft,Ha,Ht,CF] = fmsolver_decolagem2022(X,U,AircraftData,estado_do_aviao)
u = X(1);
w = X(3);
q = X(5);
teta = X(11);

% Vx = X[0]
% Vz = X[1]
% q = X[2]
% x_pos = X[3]
% z_pos = X[4]
% teta = X[5]

% Controle
dt = U(1);
de = U(2);

V = sqrt(u^2 + w^2);
%% Verificar fase de decolagem do avião
alfa_max = 20 * pi / 180;
if strcmp(estado_do_aviao, 'subida')
    % Jeito EDI de calcular o alfa
    % alfa = atan(w/u);
    % Jeito do pdf de calcular os angulos
    gama = atan(w/u);
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

% Dados gerais
geral = AircraftData{2,1};
rho = geral(1,1);
g = geral(1,2);
Sref = geral(1,3);
cref = geral(1,5);
nh = geral(1,6);
zp = geral(1,8);
xcg = geral(1,10);
zcg = geral(1,11);
zh = geral(1,12);

% Dados das Superfícies
superficies = AircraftData{3,1};
Sw = superficies(1,1);
CLaw = superficies(1,2);
CL0w = superficies(1,3);
Cmacw = superficies(1,4);
iw = superficies(1,7);
depsilondalpha = superficies(1,8);
epsilon0 = superficies(1,9);
CDw_c = [superficies(1,10) superficies(1,11) superficies(1,12)];
xac = superficies(1,13);
CLmax = superficies(1,14);
% Dados EH
Sh = superficies(2,1);
CLah = superficies(2,2);
CL0h = superficies(2,3);
Cmach = superficies(2,4);
% demin = superficies(2,5);
% demax = superficies(2,6);
ih = superficies(2,7);
Vh = superficies(2,8);
CDh_c = [superficies(2,10) superficies(2,11) superficies(2,12)];

motorInput = AircraftData{5,1};
T = max(tracao(V,dt,motorInput),0);
T = T*g;

% Derivadas da aeronave
derivadas = AircraftData{4,1};
% q
CLq = derivadas(3,1);
Cmq = derivadas(3,2);
% Profundor
CLde = derivadas(6,1);
Cmde = derivadas(6,2);

% Cálculo da pressão dinãmica
Q = 0.5*rho*V^2;

Ft = [T 0 0];
Mt = T*zp;
Ht = [0 Mt 0];

alfah = alfa*(1-depsilondalpha) - epsilon0 + ih;
CDw = CDw_c(1)*(alfa+iw)^2 + CDw_c(2)*(alfa+iw) + CDw_c(3);
CDh = CDh_c(1)*alfah^2 + CDh_c(2)*alfah + CDh_c(3);
CD = CDw + nh*Sh/Sw*CDh;
% EQUAÇÃO DE FORÇA EM X
D = Q*Sref*CD;
CL = CLaw*(alfa+iw) + CL0w + nh*Sh/Sw*(CLah*alfah + CL0h) + CLde*de;
if V ~= 0
    CL = CL + CLq*q*cref/(2*V); 
end

FS = 0.95;
if CL > FS*CLmax
    CL = CLmax;
end
L = Q*Sref*CL;

Fa = [D 0 L];

Cmcgw = (xcg-xac)*((CL0w + CLaw*(alfa+iw))*cos(alfa+iw)+CDw*sin(alfa+iw)) + zcg/cref*((CL0w + CLaw*(alfa+iw))*sin(alfa+iw)-CDw*cos(alfa+iw)) + Cmacw;
Vhz = zh*Sh/(cref*Sref);
Cmcgh = -nh*Vh*((CL0h+CLah*alfah)*cos(alfah)+CDh*sin(alfah)) -nh*Vhz*(CDh*cos(alfah)-(CL0h+CLah*alfah)*sin(alfah))+ nh*Sh/Sw*Cmach;
Cm =  Cmcgw + Cmcgh + Cmde*de;
if V ~= 0
    Cm = Cm + Cmq*q/(2*V);
end
Ma = Q*Sref*cref*Cm;

Ha = [0 Ma 0];
CF = [CD 0 CL 0 Cm 0];
end