% code for calculating constants

clc;
clear;

fun = @vasuu;

e = 0.82;
Ucond = 46;
K = 0.001;
Uloss = 11.03;
f = 0.94;

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

P0 = 7.384*10^3;
A = 67.35;
B = -7218.15;
C = -7.9939;
D = 0.00052333;
Ma = 0.028966;
Mw = 0.018016;

% Enter temperatures ---------------------------------------

T = [47.4,68.9,46.5,43.4,49.7]+273*[1,1,1,1,1];
F = zeros(4,1);
F(1) = G*((c(5)*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma)))+L*Cp*(T1-T(1))- (c(4))*Ac*(((T(4)+T(5))/2)-Tamb));
F(2) = L*Cp*(T(1)-T1)-e*(c(2))*Acond*(T(5)-T(1)-T(4)+T1)/(log((T(5)-T(1))/(T(4)-T1)));
F(3) = G*((Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma))-(c(5)*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma)))+L*Cp*(T(2)-T(3))-(c(4))*Ae*(((T(4)+T(5))/2)-Tamb));
F(4) = G*((c(5)*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma))) - e*(c(3))*a*V*(Cp*T(2)-(c(5)*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cp*T(3))+(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma)))/log((Cp*T(2)-(c(5)*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma)))/((Cp*T(3))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma)))))));


fprintf('fun1 is: %i\n: ',F(1));
fprintf('fun2 is: %i\n: ',F(2));
fprintf('fun3 is: %i\n: ',F(3));
fprintf('fun4 is: %i\n: ',F(4));

end

