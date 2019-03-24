function [P_ehd,deltaP_ehd,h_ehd, T_ehd] = rad_s_outhead(P_out_tube, h_out_tube, tube_pass, G_tube, geometry, ref)
%% ���ڼ��ܼ���
% input parameters: ������·�������ֲ�
% output parametes:  
gravity=geometry(21);
g=geometry(22);
L_ehd=geometry(27);
D_ehd=geometry(28);
A_ehd=geometry(29);
k_M_ehd=geometry(30);
direction=geometry(33); %��������������

%% 
    P_ehd(1)=P_out_tube(1); % �ܳ�1����ѹ��
    for TN=1:1:(tube_pass-1)%�ӱ��1���㵽��һ��
        mr_ehd(TN)=sum(G_tube(1:TN)); % 
        h_ehd(TN)=h_out_tube(1:TN)*G_tube(1:TN)/sum(G_tube(1:TN)); % ��������ƽ����ֵ
        G_ehd(TN)=mr_ehd(TN)/A_ehd; % ƽ����������
        Dr_ehd=D_ehd;
        viscosity_out=CoolProp.PropsSI('V','P',P_ehd(TN)*1000,'H',h_ehd(TN)*1000,ref);
        rho_r=CoolProp.PropsSI('D','P',P_ehd(TN)*1000,'H',h_ehd(TN)*1000,ref);
        T_ehd(TN)=CoolProp.PropsSI('T','P',P_ehd(TN)*1000,'H',h_ehd(TN)*1000,ref)-273.15;
        Re_ehd=(G_ehd(TN))*Dr_ehd/viscosity_out;
        rho_ehd(TN)=rho_r;
            if (Re_ehd<=2500)
                f1_ehd=16/Re_ehd;
            elseif (Re_ehd<=20000)&&(Re_ehd>2500)
                f1_ehd=0.079*Re_ehd^(-0.25);
            else
                f1_ehd=0.046*Re_ehd^(-0.2);
            end
        %���㵥���ܵ�ѹ���ͻ�����������״̬�Ȳ���
        %���ڼ���ѹ����FĦ��ѹ����D���ڼ���ѹ����G����ѹ����M��ܲ����ֲ�ѹ��
        deltaP_F_ehd(TN)=f1_ehd*4*G_ehd(TN)^2*L_ehd/(2*rho_r*Dr_ehd)*0.001; % Ħ��ѹ��
        deltaP_D_ehd(TN)=0;%G_ehd(TN)^2/rho_r*0.001;
        deltaP_G_ehd(TN)=direction*gravity*L_ehd*rho_ehd(TN)*g*0.001; % ����ѹ��
        deltaP_M_ehd(TN)=k_M_ehd*G_ehd(TN)^2/(2*rho_ehd(TN))*0.001;
        deltaP_ehd(TN)=deltaP_F_ehd(TN)+deltaP_D_ehd(TN)+deltaP_G_ehd(TN)+deltaP_M_ehd(TN); % ��ѹ��
        P_ehd(TN+1)=P_ehd(TN)-deltaP_ehd(TN);
    end
deltaP_ehd(tube_pass)=0;
h_ehd(tube_pass)=h_out_tube*G_tube/sum(G_tube); % ���������ƽ����ֵ
T_ehd(tube_pass)=CoolProp.PropsSI('T','P',P_ehd(tube_pass)*1000,'H',h_ehd(tube_pass)*1000,ref)-273.15;
end