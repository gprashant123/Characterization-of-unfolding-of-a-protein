function[yval,pu,pf,H,S,G] = myExp(init,T)

%{

This function takes a vector of parameters and the temperature as
arguments. It then calculates the value of the signal using thermodynamic
equations.

%}

R=8.314;
Sm=init(2)/init(1);
H=init(2) + init(3)*(T-init(1));
S=Sm + init(3)*log(T/init(1));
G= init(2) + init(3)*(T-init(1)) - T.*(Sm + init(3)*log(T/init(1)));
Keq= exp(-G./(R*T));
pf= 1./(1+Keq);
pu= Keq./(1+Keq);

sf= init(4)*T + init(5);
su= init(6)*T + init(7);
yval = sf.*pf + su.*pu;

end

