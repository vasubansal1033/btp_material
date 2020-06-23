% first version of code

clc;
clear;


% P0i = (P0*exp(A+(B/Ti)+(C*log(Ti))+D*Ti));
% Hi = (Cpa*Ti + (Cpw*Ti+Hvap)*((P0*exp(A+(B/Ti)+(C*log(Ti))+D*Ti))*Mw)/((P-(P0*exp(A+(B/Ti)+(C*log(Ti))+D*Ti)))*Ma));
% H1 = (Cpa*T1 + (Cpw*T1+Hvap)*((P0*exp(A+(B/T1)+(C*log(T1))+D*T1))*Mw)/((P-(P0*exp(A+(B/T1)+(C*log(T1))+D*T1)))*Ma));
% H2 = (Cpa*T(1) + (Cpw*T(1)+Hvap)*((P0*exp(A+(B/T(1))+(C*log(T(1)))+D*T(1)))*Mw)/((P-(P0*exp(A+(B/T(1))+(C*log(T(1)))+D*T(1))))*Ma));
% H3 = (Cpa*T(2) + (Cpw*T(2)+Hvap)*((P0*exp(A+(B/T(2))+(C*log(T(2)))+D*T(2)))*Mw)/((P-(P0*exp(A+(B/T(2))+(C*log(T(2)))+D*T(2))))*Ma));
% H4 = (Cpa*T(3) + (Cpw*T(3)+Hvap)*((P0*exp(A+(B/T(3))+(C*log(T(3)))+D*T(3)))*Mw)/((P-(P0*exp(A+(B/T(3))+(C*log(T(3)))+D*T(3))))*Ma));
% H5 = (Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma));
% H6 = (Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma));

fun = @vasuu;

T0 = [71,71,24,60,82]+[1,1,1,1,1]*273;
opts = optimoptions('fsolve', 'TolFun', 1E-40, 'TolX', 1E-40);
T = fsolve(fun,T0,opts);

know = T-[273,273,273,273,273]
% D_ = G*((f*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*(Mw/Ma)/(P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))))-(((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*(Mw/Ma)/(P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))))))*3600
%H6= ((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma)
%H5= ((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma)
%H4= ((P0*exp(A+(B/T(3))+(C*log(T(3)))+D*T(3)))*Mw)/((P-(P0*exp(A+(B/T(3))+(C*log(T(3)))+D*T(3))))*Ma)
%H3= ((P0*exp(A+(B/T(2))+(C*log(T(2)))+D*T(2)))*Mw)/((P-(P0*exp(A+(B/T(2))+(C*log(T(2)))+D*T(2))))*Ma)

function F = vasuu(T)

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
Uloss = 7.04;
Ulc = Uloss;
Ule = Uloss;
Ac = 2*(0.305*0.335*2+0.305*0.305);
Ae = 2*(0.305*0.335*2+0.305*0.305);
e = 0.82;
Ucond = 47.90;
K = 0.0015;
f = 0.82;

P0 = 7.384*10^3;
A = 67.35;
B = -7218.15;
C = -7.9939;
D = 0.00052333;
Ma = 0.028966;
Mw = 0.018016;


F = zeros(5,1);
F(1) = G*(f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma)))+L*Cp*(T1-T(1))- Ulc*Ac*(((T(4)+T(5))/2)-Tamb);
F(2) = L*Cp*(T(1)-T1)-e*Ucond*Acond*(T(5)-T(1)-T(4)+T1)/(log((T(5)-T(1))/(T(4)-T1)));
F(3) = G*((Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma))-f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma)))+L*Cp*(T(2)-T(3))-Ule*Ae*(((T(4)+T(5))/2)-Tamb);
F(4) = G*(f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma))) - e*K*a*V*(Cp*T(2)-f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cp*T(3))+(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma)))/log((Cp*T(2)-f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma)))/((Cp*T(3))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma))));
F(5) = L*Cp*(T(2)-T(1))-Qdot;


fprintf('fun1 is: %i\n: ',F(1));
fprintf('fun2 is: %i\n: ',F(2));
fprintf('fun3 is: %i\n: ',F(3));
fprintf('fun4 is: %i\n: ',F(4));
fprintf('fun5 is: %i\n: ',F(5));

end

