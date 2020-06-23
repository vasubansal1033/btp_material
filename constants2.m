% code for calculating constants

clc;
clear;

fun = @vasuu;

e = 0.82;
Ucond = 46.78;
K = 0.0014;
Uloss = 10.7174;
f = 0.8158;

c0 = [0.82,Ucond,K,Uloss,f];
opts = optimoptions('fsolve', 'TolFun', 1E-20, 'TolX', 1E-20);
c = fsolve(fun,c0,opts);
c
function F = vasuu(c)

% Constants
L=0.020;
G=0.040;
Cpw=4.193*10^3;
Cpa=1.009*10^3;
Hvap = 2332.20*10^3;
P = 101.325*10^3;
Acond = 3.5;
Tamb = 28+273;
T1 = 30+273;
Cp = 2.4359*10^3;
V = 0.0143;
a = 100;
% A = 2*(0.305*0.335*2+0.305*0.335);
Qdot = 1120;
% Uloss = 7.04;
% Ulc = Uloss;
% Ule = Uloss;
Ac = 2*(0.305*0.335*2+0.305*0.305);
Ae = 2*(0.305*0.335*2+0.305*0.305);
e = 0.82;
% Ucond = 47.90;
% K = 0.0015;
% f = 0.82;

Ma = 0.028966;
Mw = 0.018016;


% Enter temperatures ---------------------------------------

T = [47.1699,70.1594,48.1215,43.9680,49.2383]+273.15*[1,1,1,1,1];
% T = [45.7858,68.7753,47.0392,42.5009,47.9268]+273.15*[1,1,1,1,1];
H3 = @(T) (Cp*T(2));
H4 = @(T) (Cp*T(3));
H5 = @(T) (Cpa*T(4) + (Cpw*T(4)+Hvap)*610.78*exp((17.27*(T(4)-273))/(T(4)-273+237.3))*(Mw/Ma))/(P-(610.78*exp((17.27*(T(4)-273))/(T(4)-273+237.3))));
H6 = @(T) (Cpa*T(5) + (Cpw*T(5)+Hvap)*610.78*exp((17.27*(T(5)-273))/(T(5)-273+237.3))*(Mw/Ma))/(P-(610.78*exp((17.27*(T(5)-273))/(T(5)-273+237.3))));

F = zeros(4,1);
F(1) = G*((c(5))*H6(T)-H5(T))+L*Cp*(T1-T(1))-(c(4))*Ac*(((T(4)+T(5))/2)-Tamb);
F(2) = L*Cp*(T(1)-T1) - (e)*(c(2))*Acond*(T(5)-T(1)-T(4)+T1)/(log((T(5)-T(1))/(T(4)-T1)));
F(3) = G*(H5(T)-(c(5))*H6(T))+L*Cp*(T(2)-T(3))-(c(4))*Ae*(((T(4)+T(5))/2)-Tamb);
F(4) = G*((c(5))*H6(T)-H5(T))-(e)*(c(3))*a*V*((H3(T)-(c(5))*H6(T)-H4(T)+H5(T))/log((H3(T)-(c(5))*H6(T))/(H4(T)-H5(T))));


fprintf('fun1 is: %i\n: ',F(1));
fprintf('fun2 is: %i\n: ',F(2));
fprintf('fun3 is: %i\n: ',F(3));
fprintf('fun4 is: %i\n: ',F(4));

lolol = L*Cp*(T1-T(1)+T(2)-T(3))/(2*Ae*(((T(4)+T(5))/2)-Tamb));
fprintf('Uloss is: %i\n',lolol);

lolol2 = (L*Cp*(T(1)-T1))/((e)*Acond*(T(5)-T(1)-T(4)+T1)/(log((T(5)-T(1))/(T(4)-T1))));
fprintf('Ucond is: %i\n',lolol2);

end

