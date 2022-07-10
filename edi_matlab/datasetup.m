function [data,setup] = datasetup()
data = cell(6,1);
controle = xlsread('dinamicadados.xlsx','Controle');
tipo = controle(1,1);
V_0 = controle(2,1);
alfa_0 = controle(3,1);
beta_0 = controle(4,1);
m = controle(5,1);
h0 = controle(6,1);
gama = controle(7,1)*pi/180;
n = controle(8,1);
data{1,1} = [tipo V_0 alfa_0 beta_0 m h0 gama n];
setup = data{1,1};

geral = xlsread('dinamicadados.xlsx','Geral');
rho = geral(1,1);
g = geral(2,1);
Sref = geral(3,1);
bref = geral(4,1);
cref = geral(5,1);
nh = geral(6,1);
yp = geral(7,1);
zp = geral(8,1);
lp = geral(9,1)*pi/180;
xcg = geral(10,1);
zcg = geral(11,1);
zh = geral(12,1);
Ib = [geral(13,1) geral(13,2) geral(13,3); geral(14,1) geral(14,2) geral(14,3); geral(15,1) geral(15,2) geral(15,3)];
x_tdp = geral(16,1);
x_tdn = geral(17,1);
mi = geral(18,1);
Sd = geral(19,1);
z_inicial = geral(20,1);
data{7,1} = Ib;
data{2,1} = [rho g Sref bref cref nh yp zp lp xcg zcg zh x_tdp x_tdn mi Sd z_inicial];

superficies = xlsread('dinamicadados.xlsx','Superficies');
% Dados asa
Sw = superficies(1,1);
CLaw = superficies(2,1)*180/pi;
CL0w = superficies(3,1);
Cmacw = superficies(4,1);
damin = superficies(5,1)*pi/180;
damax = superficies(6,1)*pi/180;
iw = superficies(7,1)*pi/180;
depsilondalpha = superficies(8,1);
epsilon0 = superficies(9,1)*pi/180;
xac = superficies(10,1);
CLmax = superficies(11,1);
% Dados EH
Sh = superficies(1,3);
CLah = superficies(2,3)*180/pi;
CL0h = superficies(3,3);
Cmach = superficies(4,3);
demin = superficies(5,3)*pi/180;
demax = superficies(6,3)*pi/180;
ih = superficies(7,3)*pi/180;
Vh = superficies(8,3);
% Dados EV
drmin = superficies(5,5)*pi/180;
drmax = superficies(6,5)*pi/180;
% Dados miscelâneos
dtmin = superficies(1,7);
dtmax = superficies(2,7);
nhelice = superficies(3,7);
data{3,1} = [Sw CLaw CL0w Cmacw damin damax iw depsilondalpha epsilon0 xac CLmax;...
             Sh CLah CL0h Cmach demin demax ih Vh 0 0 0;...
             0 0 0 0 drmin drmax 0 0 0 0 0;...
             0 0 0 0 dtmin dtmax 0 0 0 0 0];
% Arrastos
arrastos = xlsread('dinamicadados.xlsx','Arrastos');
CDpw_c1 = arrastos(1,1);
CDpw_c2 = arrastos(2,1);
CDiw_c1 = arrastos(4,1)*(180/pi)^2;
CDiw_c2 = arrastos(5,1)*(180/pi);
CDiw_c3 = arrastos(6,1);
CDph_c1 = arrastos(1,4);
CDph_c2 = arrastos(2,4);
CDih_c1 = arrastos(4,4)*(180/pi)^2;
CDih_c2 = arrastos(5,4)*(180/pi);
CDih_c3 = arrastos(6,4);
CDpv_c1 = arrastos(1,7);
CDpv_c2 = arrastos(2,7);
CDpf_c1 = arrastos(1,10);
CDpf_c2 = arrastos(2,10);
data{4,1} = [
            CDpw_c1 CDpw_c2 0; ...
            CDiw_c1 CDiw_c2 CDiw_c3; ...
            CDph_c1 CDph_c2 0; ...
            CDih_c1 CDih_c2 CDih_c3; ...
            CDpv_c1 CDpv_c2 0; ...
            CDpf_c1 CDpf_c2 0
            ];

% Derivadas
derivadas = xlsread('dinamicadados.xlsx','Derivadas');
CYb = derivadas(1,1);
Clb = derivadas(2,1);
Cnb = derivadas(3,1);
CYp = derivadas(4,1);
Clp = derivadas(5,1);
Cnp = derivadas(6,1);
CLq = derivadas(7,1);
Cmq = derivadas(8,1);
CYr = derivadas(9,1);
Clr = derivadas(10,1);
Cnr = derivadas(11,1);
CYda = derivadas(12,1);
Clda = derivadas(13,1);
Cnda = derivadas(14,1);
CLde = derivadas(15,1);
Cmde = derivadas(16,1);
CYdr = derivadas(17,1);
Cldr = derivadas(18,1);
Cndr = derivadas(19,1);
data{5,1} = [CYb Clb Cnb; CYp Clp Cnp; CLq Cmq 0; CYr Clr Cnr; CYda Clda Cnda; CLde Cmde 0; CYdr Cldr Cndr];

tracao = xlsread('dinamicadados.xlsx','Tracao');
motorInput = zeros(1,22);
for i=1:size(motorInput,2)
    if i ~= 22
        motorInput(1,i) = tracao(nhelice,i);
    else
        motorInput(1,i) = tracao(nhelice,i+2);
    end
end
data{6,1} = motorInput;
end