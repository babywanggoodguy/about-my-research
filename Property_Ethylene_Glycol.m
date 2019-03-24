function [output]=Property_Ethylene_Glycol(t)
% 输入温度单位为℃
% 此为体积分数为50%的乙二醇溶液的物性参数多项式拟合
%% Pressure polynomial
% P_coe=[1.4528E-8 -1.69953E-6 1.73105E-4 -0.00568 0.09592];
% P_t=polyval(P_coe,t);
%% Density polynomial(kg/m3)
D_coe=[-0.0024 -0.34231 1081.18652];
D_t=polyval(D_coe,t);
%% Enthalpy polynomial
% h_coe=[2.33683E-6 -7.83178E-4 0.09048 0.02193 65.28345];
% h_t=polyval(h_coe,t);
%% Cp polynomial(kJ/kg*K)
Cp_coe=[-6.23E-19 0.00386 3.2034];
Cp_t=polyval(Cp_coe,t);
%% CONDUCTIVITY polynomial(W/m*K)
Lamda_coe=[-3.75E-06 8.92E-04 0.36409];
Lamda_t=polyval(Lamda_coe,t);
%% Dynamics Viscosity polynomial(mPa*s)
u_coe=[-3.91E-06 0.00122 -0.13625 6.03901];
u_t=polyval(u_coe,t);
%% Viscosity polynomial(m2/s)
v_t=u_t/D_t;
%% Prandle number polynomial
% Pr_coe=[1.89394E-8 -1.00602E-5 0.00208 -0.20748 9.87023];
% Pr_t=polyval(Pr_coe,t);
%% output parameters
output(1)=D_t;
output(2)=Cp_t;
output(3)=Lamda_t;
output(4)=v_t;

