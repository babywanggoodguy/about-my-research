function error=rad_G_tube_sovle(G_tube0, G_pass, tubes_pass, P_inpass0, h_inlet0, va, ta, z, geometry, ref)
% 该函数是基于压力平衡得到各个管路的流量分布
%% dP in inlet-header 
[P_ihd,deltaP_ihd] = rad_s_inhead(P_inpass0, h_inlet0, tubes_pass, G_tube0, geometry, ref);
%% dP and Heat transfer in tube
P_in_tube=P_ihd;
h_in_tube=h_inlet0;
for i=1:tubes_pass
    [output,output_unit] = rad_s_tube(z, h_in_tube, P_in_tube(i), G_tube0(i), va(i,:), ta(i,:), geometry, ref);
    P_out_tube(i)=output(1);
    h_out_tube(i)=output(2);
    Qt(i)=output(3);
    deltaP_tube(i)=output(4);
    dPc_intube(i)=output(5);
    dPe_outtube(i)=output(6);
    T_out_air(i,:)=output_unit(1,:);
    deltaP_a(i,:)=output_unit(2,:);
    U(i,:)=output_unit(3,:);
    k(i,:)=output_unit(4,:);
    j_a(i,:)=output_unit(5,:);
    f_a(i,:)=output_unit(6,:);
end
%% dP in outlet-header
[P_ehd,deltaP_ehd,h_ehd,T_ehd] = rad_s_outhead(P_out_tube, h_out_tube, tubes_pass, G_tube0, geometry, ref);
T_out_rad=T_ehd(tubes_pass);
dP_ref=deltaP_tube(1)+sum(deltaP_ehd);
dP_air=mean(mean(deltaP_a));
%% error calculate
error(1)=(sum(G_tube0)-G_pass)./G_pass;
for i=2:tubes_pass
    error(i)=2.*((deltaP_ihd(i-1)+deltaP_tube(i))-(deltaP_ehd(i-1)+deltaP_tube(i-1)))./((deltaP_ihd(i-1)+deltaP_tube(i))+(deltaP_ehd(i-1)+deltaP_tube(i-1)));
end
error
