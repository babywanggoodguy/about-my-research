function [P_ehd,deltaP_ehd,h_ehd, T_ehd] = rad_s_outhead(P_out_tube, h_out_tube, tube_pass, G_tube, geometry, ref)
%% 出口集管计算
% input parameters: 各个管路的流量分布
% output parametes:  
gravity=geometry(21);
g=geometry(22);
L_ehd=geometry(27);
D_ehd=geometry(28);
A_ehd=geometry(29);
k_M_ehd=geometry(30);
direction=geometry(33); %集管内流动方向

%% 
    P_ehd(1)=P_out_tube(1); % 管程1出口压力
    for TN=1:1:(tube_pass-1)%从编号1计算到另一端
        mr_ehd(TN)=sum(G_tube(1:TN)); % 
        h_ehd(TN)=h_out_tube(1:TN)*G_tube(1:TN)/sum(G_tube(1:TN)); % 质量流量平均焓值
        G_ehd(TN)=mr_ehd(TN)/A_ehd; % 平均质量流量
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
        %计算单根管的压降和换热量及出口状态等参数
        %进口集管压降：F摩擦压降，D进口减速压降，G重力压降，M扁管插入后局部压降
        deltaP_F_ehd(TN)=f1_ehd*4*G_ehd(TN)^2*L_ehd/(2*rho_r*Dr_ehd)*0.001; % 摩擦压降
        deltaP_D_ehd(TN)=0;%G_ehd(TN)^2/rho_r*0.001;
        deltaP_G_ehd(TN)=direction*gravity*L_ehd*rho_ehd(TN)*g*0.001; % 重力压降
        deltaP_M_ehd(TN)=k_M_ehd*G_ehd(TN)^2/(2*rho_ehd(TN))*0.001;
        deltaP_ehd(TN)=deltaP_F_ehd(TN)+deltaP_D_ehd(TN)+deltaP_G_ehd(TN)+deltaP_M_ehd(TN); % 总压降
        P_ehd(TN+1)=P_ehd(TN)-deltaP_ehd(TN);
    end
deltaP_ehd(tube_pass)=0;
h_ehd(tube_pass)=h_out_tube*G_tube/sum(G_tube); % 求出口质量平均焓值
T_ehd(tube_pass)=CoolProp.PropsSI('T','P',P_ehd(tube_pass)*1000,'H',h_ehd(tube_pass)*1000,ref)-273.15;
end