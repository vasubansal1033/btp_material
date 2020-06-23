% evaluating equations of the two versions for set of Ts

clc;
clear;

% Constants
L=0.020;
G=0.040;
Cpw=4.193*10^3;
Cpa=1.009*10^3;
Hvap = 2332.20*10^3;
P = 101.325*10^3;
Acond = 3.5;
Tamb = 30.5+273;
T1 = 27.9+273;
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

T=[47.4,68.9,46.5,43.4,49.7]+273*[1,1,1,1,1];

F = zeros(5,1);
F(1) = G*(f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma)))+L*Cp*(T1-T(1))- Ulc*Ac*(((T(4)+T(5))/2)-Tamb);
F(2) = L*Cp*(T(1)-T1)-e*Ucond*Acond*(T(5)-T(1)-T(4)+T1)/(log((T(5)-T(1))/(T(4)-T1)));
F(3) = G*((Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma))-f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma)))+L*Cp*(T(2)-T(3))-Ule*Ae*(((T(4)+T(5))/2)-Tamb);
F(4) = G*(f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma))) - e*K*a*V*(Cp*T(2)-f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cp*T(3))+(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma)))/log((Cp*T(2)-f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma)))/((Cp*T(3))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma))));
F(5) = L*Cp*(T(2)-T(1))-Qdot;

fprintf('f1 is: %i\n: ',F(1));
fprintf('f2 is: %i\n: ',F(2));
fprintf('f3 is: %i\n: ',F(3));
fprintf('f4 is: %i\n: ',F(4));
fprintf('f5 is: %i\n: ',F(5));

H3 = @(T) (Cp*T(2));
H4 = @(T) (Cp*T(3));
H5 = @(T) (Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*(Mw/Ma))/(P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))));
H6 = @(T) (Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*(Mw/Ma))/(P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))));


fun = zeros(5,1);
fun(1) = G*(f*H6(T)-H5(T))+L*Cp*(T1-T(1))-Ulc*Ac*(((T(4)+T(5))/2)-Tamb);
fun(2) = L*Cp*(T(1)-T1) - e*Ucond*Acond*(T(5)-T(1)-T(4)+T1)/(log((T(5)-T(1))/(T(4)-T1)));
fun(3) = G*(H5(T)-f*H6(T))+L*Cp*(T(2)-T(3))-Ule*Ae*(((T(4)+T(5))/2)-Tamb);
fun(4) = G*(f*H6(T)-H5(T))-e*K*a*V*((H3(T)-f*H6(T)-H4(T)+H5(T))/log((H3(T)-f*H6(T))/(H4(T)-H5(T))));
fun(5) = -Qdot+L*Cp*(T(2)-T(1));

fprintf('fun1 is: %i\n: ',fun(1));
fprintf('fun2 is: %i\n: ',fun(2));
fprintf('fun3 is: %i\n: ',fun(3));
fprintf('fun4 is: %i\n: ',fun(4));
fprintf('fun5 is: %i\n: ',fun(5));

% F(3) = G*(f*H6r-H5r) - e*K*a*V*(H3r-f*H6r-H4r+H5r)/log((H3r-f*H6r)/(H4r-H5r));
% F(4) = G*(f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma))) - e*K*a*V*((Cp*T(2))-f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma))-(Cpa*T(3) + (Cpw*T(3)+Hvap)*((P0*exp(A+(B/T(3))+(C*log(T(3)))+D*T(3)))*Mw)/((P-(P0*exp(A+(B/T(3))+(C*log(T(3)))+D*T(3))))*Ma))+(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma)))/log(((Cp*T(2))-f*(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))*Ma)))/((Cpa*T(3) + (Cpw*T(3)+Hvap)*((P0*exp(A+(B/T(3))+(C*log(T(3)))+D*T(3)))*Mw)/((P-(P0*exp(A+(B/T(3))+(C*log(T(3)))+D*T(3))))*Ma))-(Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4))))*Ma))));
% fprintf('new F4 is: %i\n: ',F(4));