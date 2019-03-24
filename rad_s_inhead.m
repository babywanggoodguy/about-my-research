function [P_ihd,deltaP_ihd] = rad_s_inhead(P_in_pass, h_in_pass, tube_pass, G_tube, geometry, ref )
%进口集管计算
% 输入参数为：质量流量分布、散热器几何结构
% 输出参数为：通过集管后的压力及过程压降。
gravity=geometry(21);
g=geometry(22);
L_ihd=geometry(23);
D_ihd=geometry(24);
A_ihd=geometry(25);
k_M_ihd=geometry(26);
direction=geometry(33);
P_ihd(1)=P_in_pass;
    for TN=1:(tube_pass-1) 
        mr_ihd(TN)=sum(G_tube(TN:tube_pass)); 
        G_ihd(TN)=mr_ihd(TN)/A_ihd;
        Dr_ihd=D_ihd;
        viscosity_in=CoolProp.PropsSI('V','P',P_ihd(TN)*1000,'H',h_in_pass*1000,ref);
        rho_r=CoolProp.PropsSI('D','P',P_ihd(TN)*1000,'H',h_in_pass*1000,ref);
        Re_ihd=(G_ihd(TN))*Dr_ihd/viscosity_in;
        rho_ihd(TN)=rho_r;
            if (Re_ihd<=2500)
                f1_ihd=16/Re_ihd;
            elseif (Re_ihd<=20000)&&(Re_ihd>2500)
                f1_ihd=0.079*Re_ihd^(-0.25);
            else
                f1_ihd=0.046*Re_ihd^(-0.2);
            end
        %进口集管压降：F摩擦压降，D进口减速压降，G重力压降，M扁管插入后局部压降 (Yin J M,2002)
        deltaP_F_ihd(TN)=f1_ihd*4*G_ihd(TN)^2*L_ihd/(2*rho_r*Dr_ihd)*0.001;
        deltaP_D_ihd(TN)=0;%G_ihd(m)^2/rho_r*0.001;
        deltaP_G_ihd(TN)=direction*gravity*L_ihd*rho_ihd(TN)*g*0.001;
        deltaP_M_ihd(TN)=k_M_ihd*G_ihd(TN)^2/(2*rho_ihd(TN))*0.001;
        deltaP_ihd(TN)=deltaP_F_ihd(TN)+deltaP_D_ihd(TN)+deltaP_G_ihd(TN)+deltaP_M_ihd(TN);
        P_ihd(TN+1)=P_ihd(TN)-deltaP_ihd(TN);
    end
deltaP_ihd(tube_pass)=0;
end