% second version of second version of code

clear;
clc;

calc_temp;
function calc_temp

T = @(T2,T3,T4,T5,T6) [T2 T3 T4 T5 T6]

% Constants
L=0.20;
G=0.040;
Cpw=4.193*10^3;
Cpa=1.009*10^3;
Hvap = 2332.20*10^3;
P = 101.325*10^3;
Acond = 3.5;
Tamb = 30.5+273;
T1 = 27.9+273;
Cp = 2.4359*10^3;
% Cp = 3.12*10^3;
V = 0.0143;
a = 100;
% A = 2*(0.305*0.335*2+0.305*0.335);
Qdot = 1000;
Uloss = 7.04;
Ulc = Uloss;
Ule = Uloss;
Ac = 2*(0.305*0.335*2+0.305*0.305);
Ae = 2*(0.305*0.335*2+0.305*0.305);
e = 0.82;
Ucond = 47.90;
K = 0.0015;
f = 0.82;


Ma = 0.028966;
Mw = 0.018016;

% Constants obtained through inverse problem
% e = 0.82;
% Ucond = 46.78;
% K = 0.0014;
% Uloss = 10.7174;
% Ulc = Uloss;
% Ule = Uloss;
% f = 0.8158;


%------------------------

H3 = @(T) (Cp*T(2));
H4 = @(T) (Cp*T(3));


% H3 = @(T) (Cpa*T(2) + (Cpw*T(2)+Hvap)*610.78*exp((17.27*(T(2)-273))/(T(2)-273+237.3))*(Mw/Ma))/(P-(610.78*exp((17.27*(T(2)-273))/(T(2)-273+237.3))));
% H4 = @(T) (Cpa*T(3) + (Cpw*T(3)+Hvap)*610.78*exp((17.27*(T(3)-273))/(T(3)-273+237.3))*(Mw/Ma))/(P-(610.78*exp((17.27*(T(3)-273))/(T(3)-273+237.3))));

H5 = @(T) (Cpa*T(4) + (Cpw*T(4)+Hvap)*610.78*exp((17.27*(T(4)-273))/(T(4)-273+237.3))*(Mw/Ma))/(P-(610.78*exp((17.27*(T(4)-273))/(T(4)-273+237.3))));
H6 = @(T) (Cpa*T(5) + (Cpw*T(5)+Hvap)*610.78*exp((17.27*(T(5)-273))/(T(5)-273+237.3))*(Mw/Ma))/(P-(610.78*exp((17.27*(T(5)-273))/(T(5)-273+237.3))));

% Solution
opts = optimoptions('fsolve', 'TolFun', 1E-40, 'TolX', 1E-40);
T0 = [47.16,70.16,48.1215,43.968,49.2383]+273.15*[1,1,1,1,1];
Temp = fsolve(@CalcTemps, T0,opts);



function fun = CalcTemps(T)
    
    fun(1) = G*(f*H6(T)-H5(T))+L*Cp*(T1-T(1))-Ulc*Ac*(((T(4)+T(5))/2)-Tamb);
    fun(2) = L*Cp*(T(1)-T1) - e*Ucond*Acond*(T(5)-T(1)-T(4)+T1)/(log((T(5)-T(1))/(T(4)-T1)));
    fun(3) = G*(H5(T)-f*H6(T))+L*Cp*(T(2)-T(3))-Ule*Ae*(((T(4)+T(5))/2)-Tamb);
    fun(4) = G*(f*H6(T)-H5(T))-e*K*a*V*((H3(T)-f*H6(T)-H4(T)+H5(T))/log((H3(T)-f*H6(T))/(H4(T)-H5(T))));
    fun(5) = -Qdot+L*Cp*(T(2)-T(1));


W5r = 610.78*exp((17.27*(T(4)-273))/(T(4)-273+237.3))*(Mw/Ma)/(P-(610.78*exp((17.27*(T(4)-273))/(T(4)-273+237.3))));
W6r = 610.78*exp((17.27*(T(5)-273))/(T(5)-273+237.3))*(Mw/Ma)/(P-(610.78*exp((17.27*(T(5)-273))/(T(5)-273+237.3))));
D_ = G*(f*W6r-W5r)*3600;

fprintf('Production rate is: %i\n', D_);
fprintf('fun1 is: %i\n ',fun(1));
fprintf('fun2 is: %i\n ',fun(2));
fprintf('fun3 is: %i\n ',fun(3));
fprintf('fun4 is: %i\n ',fun(4));
fprintf('fun5 is: %i\n ',fun(5));


end

Temp - [1,1,1,1,1]*273
end
