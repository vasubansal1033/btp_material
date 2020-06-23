% Constants
L=0.02;
G=0.04;
Cpw=4.193*10^3;
Cpa=1.009*10^3;
Hvap = 2332.20*10^3;
Patm = 101.325*10^3;
Acond = 3.5;
Tamb = 28+273;
T1 = 30+273;
Cp = 3.12*10^3;
V = 0.0143;
a = 300;
A = 2*(0.305*0.335*2+0.305*0.335);
Qdot = 1120;
Uloss = 7.04;
Ulc = Uloss;
Ule = Uloss;
Ac = A;
Ae =A;
e = 0.82;
Ucond = 47.90;
K = 0.0053;
f = 0.82;

P0 = 7.384*10^3;
A = 67.35;
B = -7218.15;
C = -7.9939;
D = 0.00052333;
Ma = 0.028966;
Mw = 0.018016;
% P0i = (P0*exp(A+(B/Ti)+(C*ln(Ti))+D*Ti));
% Hi = (Cpa*Ti + (Cpw*Ti+Hvap)*((P0*exp(A+(B/Ti)+(C*ln(Ti))+D*Ti))*Mw)/((P-(P0*exp(A+(B/Ti)+(C*ln(Ti))+D*Ti)))*Ma));
% H1 = (Cpa*T1 + (Cpw*T1+Hvap)*((P0*exp(A+(B/T1)+(C*ln(T1))+D*T1))*Mw)/((P-(P0*exp(A+(B/T1)+(C*ln(T1))+D*T1)))*Ma));
% H2 = (Cpa*T(2) + (Cpw*T(2)+Hvap)*((P0*exp(A+(B/T(2))+(C*ln(T(2)))+D*T(2)))*Mw)/((P-(P0*exp(A+(B/T(2))+(C*ln(T(2)))+D*T(2))))*Ma));
% H3 = (Cpa*T(3) + (Cpw*T(3)+Hvap)*((P0*exp(A+(B/T(3))+(C*ln(T(3)))+D*T(3)))*Mw)/((P-(P0*exp(A+(B/T(3))+(C*ln(T(3)))+D*T(3))))*Ma));
% H4 = (Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*ln(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*ln(T(4)))+D*T(4))))*Ma));
% H5 = (Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5))))*Ma));
% H6 = (Cpa*T(6) + (Cpw*T(6)+Hvap)*((P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6)))*Mw)/((P-(P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6))))*Ma));


function F = vasu(T)
F(1) = G*(f*(Cpa*T(6) + (Cpw*T(6)+Hvap)*((P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6)))*Mw)/((P-(P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6))))*Ma))-(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5))))*Ma))) + L*Cp(T1-T(1))- Ulc*Ac*(((T(4)+T(5))/2)-Tamb);
F(2) = L*Cp*(T(2)-T1)-e*Ucond*Acond*(T(5)-T(1)-T(4)+T1)/ln((T(5)-T(1))/(T(4)-T1));
F(3) = G*((Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5))))*Ma))-f*(Cpa*T(6) + (Cpw*T(6)+Hvap)*((P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6)))*Mw)/((P-(P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6))))*Ma)))+L*Cp*(T(2)-T(3))-Ule*Ae*(((T(4)+T(5))/2)-Tamb);
F(4) = G*(f*(Cpa*T(6) + (Cpw*T(6)+Hvap)*((P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6)))*Mw)/((P-(P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6))))*Ma))-(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5))))*Ma)))- e*K*a*V*((Cpa*T(3) + (Cpw*T(3)+Hvap)*((P0*exp(A+(B/T(3))+(C*ln(T(3)))+D*T(3)))*Mw)/((P-(P0*exp(A+(B/T(3))+(C*ln(T(3)))+D*T(3))))*Ma))-f*(Cpa*T(6) + (Cpw*T(6)+Hvap)*((P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6)))*Mw)/((P-(P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6))))*Ma))-((Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*ln(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*ln(T(4)))+D*T(4))))*Ma))-(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5))))*Ma))))/(ln((Cpa*T(3) + (Cpw*T(3)+Hvap)*((P0*exp(A+(B/T(3))+(C*ln(T(3)))+D*T(3)))*Mw)/((P-(P0*exp(A+(B/T(3))+(C*ln(T(3)))+D*T(3))))*Ma))-f*(Cpa*T(6) + (Cpw*T(6)+Hvap)*((P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6)))*Mw)/((P-(P0*exp(A+(B/T(6))+(C*ln(T(6)))+D*T(6))))*Ma)))/((Cpa*T(4) + (Cpw*T(4)+Hvap)*((P0*exp(A+(B/T(4))+(C*ln(T(4)))+D*T(4)))*Mw)/((P-(P0*exp(A+(B/T(4))+(C*ln(T(4)))+D*T(4))))*Ma))-(Cpa*T(5) + (Cpw*T(5)+Hvap)*((P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5)))*Mw)/((P-(P0*exp(A+(B/T(5))+(C*ln(T(5)))+D*T(5))))*Ma))));
F(5) = L*Cp*(T(2)-T(1));

fun = @vasu;
T0 = [0,0,0,0,0];
T = fsolve(fun,T0);

end