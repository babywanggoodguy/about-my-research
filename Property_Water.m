function [output]=Property_Water(t)
% 输入温度单位为℃
% 有效使用温度范围为40~140℃
%% Pressure polynomial
P_coe=[1.4528E-8 -1.69953E-6 1.73105E-4 -0.00568 0.09592];
P_t=polyval(P_coe,t)*10^5;
%% Density polynomial
D_coe=[-2.91375E-8 1.45299E-5 -0.0049 -0.05688 1001.49068];
D_t=polyval(D_coe,t);
%% Enthalpy polynomial
h_coe=[2.33683E-6 -7.83178E-4 0.09048 0.02193 65.28345];
h_t=polyval(h_coe,t);
%% Cp polynomial
Cp_coe=[1.1655E-9 -4.45221E-7 6.98718E-5 -0.00396 4.24614];
Cp_t=polyval(Cp_coe,t);
%% CONDUCTIVITY polynomial
Lamda_coe=[2.91375E-9 2.07848E-6 -0.00149 0.25204 55.66364];
Lamda_t=polyval(Lamda_coe,t)/100;
%% Viscosity polynomial
v_coe=[1.48601E-9 -9.392E-7 2.22194E-4 -0.02491 1.35562];
v_t=polyval(v_coe,t)*10^-6;
%% Prandle number polynomial
Pr_coe=[1.89394E-8 -1.00602E-5 0.00208 -0.20748 9.87023];
Pr_t=polyval(Pr_coe,t);
%% output parameters
output(1)=D_t;
output(2)=Cp_t;
output(3)=Lamda_t;
output(4)=v_t;
output(5)=P_t;
output(6)=h_t;
output(7)=Pr_t;
