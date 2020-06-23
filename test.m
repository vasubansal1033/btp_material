clear;
clc;

calc_temp
function calc_temp

T = @(T2,T3,T4,T5,T6) [T2 T3 T4 T5 T6]

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
Cp = 3.12*10^3;
V = 0.0143;
a = 300;
% A = 2*(0.305*0.335*2+0.305*0.335);
Qdot = 0;
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

%------------------------


% Solution
opts = optimoptions('fsolve', 'TolFun', 1E-30, 'TolX', 1E-30);
T0 = [28.9,31.5,30.8,30,52.89]+273*[1,1,1,1,1];
Temp = fsolve(@CalcTemps, T0,opts);

function fun = CalcTemps(T)
    H3 = @(T) (Cpa*T(2) + (Cpw*T(2)+Hvap)*((P0*exp(A+(B/T(2))+(C*log(T(2)))+D*T(2)))*(Mw/Ma))/(P-(P0*exp(A+(B/T(2))+(C*log(T(2)))+D*T(2)))));
    H4 = @(T) (Cpa*T(3) + (Cpw*T(3)+Hvap)*((P0*exp(A+(B/T(3))+(C*log(T(3)))+D*T(3)))*(Mw/Ma))/(P-(P0*exp(A+(B/T(3))+(C*log(T(3)))+D*T(3)))));
    H5 = @(T) (Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*(Mw/Ma))/(P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))));
    H6 = @(T) (Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*(Mw/Ma))/(P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))));

    fun(1) = G*(f*H6(T)-H5(T))+L*Cp*(T1-T(1))-Ulc*Ac*(((T(4)+T(5))/2)-Tamb);
    fun(2) = L*Cp*(T(1)-T1) - e*Ucond*Acond*(T(5)-T(1)-T(4)+T1)/(log((T(5)-T(1))/(T(4)-T1)));
    fun(3) = G*(H5(T)-f*H6(T))+L*Cp*(T(2)-T(3))-Ule*Ae*(((T(4)+T(5))/2)-Tamb);
    fun(4) = G*(f*H6(T)-H5(T))-e*K*a*V*((H3(T)-f*H6(T)-H4(T)+H5(T))/log((H3(T)-f*H6(T))/(H4(T)-H5(T))));
    fun(5) = -Qdot+L*Cp*(T(2)-T(1));
    
D_ = G*((f*((P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5)))*(Mw/Ma)/(P-(P0*exp(A+(B/T(5))+(C*log(T(5)))+D*T(5))))))-(((P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))*(Mw/Ma)/(P-(P0*exp(A+(B/T(4))+(C*log(T(4)))+D*T(4)))))))*3600

end

Temp - [1,1,1,1,1]*273
end
