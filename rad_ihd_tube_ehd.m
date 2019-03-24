function [output_pass,output_array] = rad_ihd_tube_ehd(mr_ihd, G_tube, mr_ehd, P_ihd, T_inlet, T_ehd0, va, ta, z, geometry)
% 输入单根管子的质量流量
gravity=geometry(21);
g=geometry(22);
L_ehd=geometry(27);
D_ehd=geometry(28);
A_ehd=geometry(29);
k_M_ehd=geometry(30);
L_ihd=geometry(23);
D_ihd=geometry(24);
A_ihd=geometry(25);
k_M_ihd=geometry(26);
direction=geometry(33); % 集管内流动方向，-1表示向下，1表示向上，0表示水平；

%% inhead dP         
G_ihd=mr_ihd/A_ihd;
Dr_ihd=D_ihd; 
prop_T=Property_Ethylene_Glycol(T_inlet); %物性参数
rho_r_ihd=prop_T(1);
viscosity_in=prop_T(4);
Re_ihd=G_ihd*Dr_ihd/viscosity_in;
if (Re_ihd<=2500)
    f1_ihd=16/Re_ihd;
elseif (Re_ihd<=20000)&&(Re_ihd>2500)
    f1_ihd=0.079*Re_ihd^(-0.25);
else
    f1_ihd=0.046*Re_ihd^(-0.2); 
end

%进口集管压降：F摩擦压降，D进口减速压降，G重力压降，M扁管插入后局部压降
deltaP_F_ihd=f1_ihd*4*G_ihd^2*L_ihd/(2*rho_r_ihd*Dr_ihd)*0.001;
deltaP_D_ihd=0;%G_ihd(m)^2/rho_r*0.001;
deltaP_G_ihd=direction.*gravity*L_ihd*rho_r_ihd*g*0.001;
% deltaP_G_ihd=0; 
deltaP_M_ihd=k_M_ihd*G_ihd^2/(2*rho_r_ihd)*0.001;
deltaP_ihd=deltaP_F_ihd+deltaP_D_ihd+deltaP_G_ihd+deltaP_M_ihd;

%% dP in tube
P_in_tube=P_ihd-deltaP_ihd; 
T_in_tube=T_inlet;
[output_tube,output_unit] = rad_s_tube(z,T_in_tube,P_in_tube,G_tube,va,ta,geometry);
P_out_tube=output_tube(1);
T_out_tube=output_tube(2); % (℃)
Qt=output_tube(3); 
deltaP_tube=output_tube(4);
% dPc_intube=output_tube(5);
% dPe_outtube=output_tube(6);
output_array=output_unit;
% T_out_air=output_unit(1,:);
% deltaP_a=output_unit(2,:);
% U=output_unit(3,:);
% NTU=output_unit(4,:);
% h_a=output_unit(5,:);
% f_a=output_unit(6,:);
% DeltaP_ref=output_unit(7,:);

%% dP in outhead
T_ehd=(T_out_tube*G_tube+T_ehd0*(mr_ehd-G_tube))/mr_ehd; % 质量流量平均Temperature 
G_ehd=mr_ehd/A_ehd; % 平均质量流量
prop_T=Property_Ethylene_Glycol(T_ehd); %物性参数
rho_r_ohd=prop_T(1);
vis_out=prop_T(4);
Re_ehd=(G_ehd)*D_ehd/vis_out;
if (Re_ehd<=2500)
    f1_ehd=16/Re_ehd;
elseif (Re_ehd<=20000)&&(Re_ehd>2500)
    f1_ehd=0.079*Re_ehd^(-0.25);
else
    f1_ehd=0.046*Re_ehd^(-0.2);
end
%计算单根管的压降和换热量及出口状态等参数
%进口集管压降：F摩擦压降，D进口减速压降，G重力压降，M扁管插入后局部压降
deltaP_F_ehd=f1_ehd*4*G_ehd^2*L_ehd/(2*rho_r_ohd*D_ehd)*0.001; % 摩擦压降
deltaP_D_ehd=0;%G_ehd^2/rho_r*0.001;
deltaP_G_ehd=direction.*gravity*L_ehd*rho_r_ohd*g*0.001; % 重力压降
% deltaP_G_ehd=0; 
deltaP_M_ehd=k_M_ehd*G_ehd^2/(2*rho_r_ohd)*0.001;
deltaP_ehd=deltaP_F_ehd+deltaP_D_ehd+deltaP_G_ehd+deltaP_M_ehd; % 总压降
P_ehd=P_out_tube-deltaP_ehd;

%% output parameters
output_pass(1)=deltaP_ihd;
output_pass(2)=deltaP_tube;
output_pass(3)=deltaP_ehd;
output_pass(4)=P_in_tube;
output_pass(5)=P_out_tube;
output_pass(6)=Qt;
output_pass(7)=T_ehd;
output_pass(8)=P_ehd;
% output_array(1,:)=T_out_air;
% output_array(2,:)=deltaP_a;
% output_array(3,:)=NTU;
% output_array(4,:)=U;
end