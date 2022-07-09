function T = tracao(V,dt,motorInput)
n = dt;
v = V;

p00 = motorInput(1);
p10 = motorInput(2);
p01 = motorInput(3);
p20 = motorInput(4);
p11 = motorInput(5);
p02 = motorInput(6);
p30 = motorInput(7);
p21 = motorInput(8);
p12 = motorInput(9);
p03 = motorInput(10);
p40 = motorInput(11);
p31 = motorInput(12);
p22 = motorInput(13);
p13 = motorInput(14);
p04 = motorInput(15);
p50 = motorInput(16);
p41 = motorInput(17);
p32 = motorInput(18);
p23 = motorInput(19);
p14 = motorInput(20);
p05 = motorInput(21);
T0 = motorInput(22);

T = p00 + p10*n + p01*v + p20*n^2 + p11*n*v + p02*v^2 + p30*n^3 + p21*n^2*v  + p12*n*v^2 + p03*v^3 + p40*n^4 ...
    + p31*n^3*v + p22*n^2*v^2 + p13*n*v^3 + p04*v^4 + p50*n^5 + p41*n^4*v + p32*n^3*v^2 + p23*n^2*v^3 +...
    p14*n*v^4 + p05*v^5;
end