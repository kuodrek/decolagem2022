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
xac = superficies(1,10);
CLmax = superficies(1,11);
% Dados EH
Sh = superficies(2,1);
CLah = superficies(2,2);
CL0h = superficies(2,3);
Cmach = superficies(2,4);
% demin = superficies(2,5);
% demax = superficies(2,6);
ih = superficies(2,7);
Vh = superficies(2,8);

% Dados de arrasto
arrasto = AircraftData{4,1};
% Asa
CDpw_c = [arrasto(1,1) arrasto(1,2)];
CDiw_c = [superficies(2,1) superficies(2,2) superficies(2,3)];
% EH
CDph_c = [arrasto(3,1) arrasto(3,2)];
CDih_c = [arrasto(4,1) arrasto(4,2) arrasto(4,3)];
% EV
CDpv_c = [arrasto(5,1) arrasto(5,2)];
% Fuselagem
CDpf_c = [arrasto(6,1) arrasto(6,2)];

motorInput = AircraftData{6,1};
T = max(tracao(V,dt,motorInput),0);
T = T*g;

% Derivadas da aeronave
derivadas = AircraftData{5,1};
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
CDiw = CDiw_c(1)*(alfa+iw)^2 + CDiw_c(2)*(alfa+iw) + CDiw_c(3);
CDpw = CDpw_c(1)*V + CDpw_c(2);
CDw = CDiw + CDpw;
CDih = CDih_c(1)*alfah^2 + CDih_c(2)*alfah + CDih_c(3);
CDph = CDph_c(1)*(V*nh) + CDph_c(2);
CDh = CDih + CDph;
CDpv = CDpv_c(1)*(V*nh) + CDpv_c(2);
CDpf = CDpf_c(1)*V + CDpf_c(2);
CD = CDw + nh*Sh/Sw*CDh + CDpv + CDpf;
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