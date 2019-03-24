function [output,output_unit] = rad_s_tube(z,T_in_tube,P_in_tube,G_t,va,ta,geometry)
% ����ı�ܼ���,����������ǵ�����ܵ�
% ���룺������ѹ������ֵ�������������ٶȡ��¶�
% �������ܳ���ѹ������ֵ����ܻ�������ˮ��ѹ��������������¶�
%% Ԥ�ȷ����ڴ棬����߼����ٶ�
IT=ones(z,1);
OT=ones(z,1);
IP=ones(z,1);
OP=ones(z,1);
Q=ones(z,1);
DeltaP_ref=ones(z,1);
h_a=ones(z,1);
AT_out=ones(z,1);
DeltaP_G=ones(z,1);
Nu_t=ones(z,1);
k=ones(z,1);
U=ones(z,1);
NTU=ones(z,1);
y_min=ones(z,1);
epsilon=ones(z,1);
deltaP_a=ones(z,1);
DeltaP_fr=ones(z,1);
f_a=ones(z,1);
output_unit=ones(z,17);

A_tr_rad=geometry(17);
sigma=geometry(32);
%% ��������ת��ѹ��

%% ��ܽ���ͻ��ѹ��
    prop_T_in=Property_Ethylene_Glycol(T_in_tube); %���Բ���
    rho_r_in=prop_T_in(1);
%     Cp_in=prop_T_in(2);
%     lambda_t_in=prop_T_in(3);
%     vis_in=prop_T_in(4);
%     Pr_in=vis_in.*rho_r_in.*Cp_in/lambda_t_in;
    rho_dpc=rho_r_in;
    Kc=0.55; % ͻ��ͻ��ѹ��ϵ���ɲ�Yin,2002
    dPc_intube=(abs(G_t)./A_tr_rad)^2./rho_dpc./2.*(1-sigma^2+Kc).*0.001; % �����ܲ�����ѹ��
    P_in_tube=P_in_tube-dPc_intube;
%% ��ܻ��ȼ�ѹ��
IT(1)=T_in_tube;
IP(1)=P_in_tube;
OP(1)=IP(1);
% Oh(1)=Ih(1).*(0.999999);
for i=1:1:z    
    volume_result=rad_s_volume(IT(i),IP(i),G_t,va(i),ta(i),geometry); % ������һ��Ԫ�ļ���   
    output_unit(i,:)=volume_result;
    Q(i)=volume_result(1); % heat transfer 
%     DeltaP_ref(i)=volume_result(2); % �ܳ�ѹ��=����ѹ��+Ħ��ѹ��
    OP(i)=volume_result(3);
%     h_a(i)=volume_result(4);
%     AT_out(i)=volume_result(5);
%     DeltaP_G(i)=volume_result(6); % ����ѹ��
%     Nu_t(i)=volume_result(7);
    OT(i)=volume_result(8);
%     k(i)=volume_result(10);
%     U(i)=volume_result(11);
%     NTU(i)=volume_result(12);
%     y_min(i)=volume_result(13);
%     epsilon(i)=volume_result(14);
%     deltaP_a(i)=volume_result(15); % ������ѹ�� Pa
%     DeltaP_fr(i)=volume_result(16); % Ħ��ѹ��
%     f_a(i)=volume_result(17);
%     err_unit(i)=volume_result(18);
    if i~=z
        IP(i+1)=OP(i);
        IT(i+1)=OT(i);
    end
end
Qt=sum(Q);
T_out_tube=OT(z);
P_out_tube=OP(z);
%% ��ܳ���ͻ��ѹ��
    prop_T_out=Property_Ethylene_Glycol(T_in_tube); %���Բ���
    rho_dpc=prop_T_out(1);
    Ke=0.99;
    dPe_outtube=(abs(G_t)./A_tr_rad)^2./rho_dpc./2.*(Ke-1+sigma^2).*0.001;
    P_out_tube=P_out_tube-dPe_outtube;
    deltaP_tube=P_in_tube-P_out_tube+dPc_intube;
%% ��������ת��ѹ��

%% outputs parameters    
output(1)=P_out_tube;
output(2)=T_out_tube; % (��)
output(3)=Qt; % heat transfer in liquid side (kW) 
output(4)=deltaP_tube; 
output(5)=dPc_intube;
output(6)=dPe_outtube;
% output_unit(1,:)=AT_out;
% output_unit(2,:)=deltaP_a; % Pressure drop in air side (Pa)
% output_unit(3,:)=U;
% output_unit(4,:)=NTU;
% output_unit(5,:)=h_a;
% output_unit(6,:)=f_a;
% output_unit(7,:)=DeltaP_ref;
end