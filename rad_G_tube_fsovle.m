function error=rad_G_tube_fsovle(G_tube0, G_pass, tube_pass, P_inpass0, T_inlet0, va, ta, z, geometry)
% 该函数是基于压力平衡得到各个管路的流量分布
% 采用二分法迭代计算真实流量值
% 调用rad_ihd_tube_ehd单管程计算子函数
%% 预先分配内存，以提高计算速度
mr_ihd=ones(tube_pass,1);
mr_ehd=ones(tube_pass,1);
P_ihd=ones(tube_pass,1);
deltaP_ihd=ones(tube_pass,1);
deltaP_tube=ones(tube_pass,1);
deltaP_ehd=ones(tube_pass,1);
P_in_tube=ones(tube_pass,1);
P_out_tube=ones(tube_pass,1);
P_ehd=ones(tube_pass,1);
T_ehd=ones(tube_pass,1);
Q_t=ones(tube_pass,1);
h_ehd=ones(tube_pass,1);
%% main calculation codes
error=ones(1,tube_pass); % 误差初始赋值
mr_ihd(1)=G_pass;
mr_ehd(1)=G_tube0(1);
P_ihd(1)=P_inpass0;
T_ehd0=0; % 初始管子出口 Assumed Temperature
[output_pass,output_array] = ...
    rad_ihd_tube_ehd(mr_ihd(1), G_tube0(1), mr_ehd(1), P_ihd(1), T_inlet0, T_ehd0, va(1,:), ta(1,:), z, geometry);
deltaP_ihd(1)=output_pass(1);
deltaP_tube(1)=output_pass(2);
deltaP_ehd(1)=output_pass(3);
P_in_tube(1)=output_pass(4);
P_out_tube(1)=output_pass(5);
Q_t(1)=output_pass(6);
T_ehd(1)=output_pass(7);
P_ehd(1)=output_pass(8);
T_out_air(1,:)=output_array(1,:);
deltaP_a(1,:)=output_array(2,:);
for i=2:tube_pass
    mr_ihd(i)=mr_ihd(i-1)-G_tube0(i);
    mr_ehd(i)=mr_ehd(i-1)+G_tube0(i);
    P_ihd(i)=P_in_tube(i-1)-deltaP_ihd(i-1);
    [output_pass,output_array] = rad_ihd_tube_ehd(mr_ihd(i), G_tube0(i), mr_ehd(i), P_ihd(i), T_inlet0, h_ehd(i-1), va(i,:), ta(i,:), z, geometry);
    deltaP_ihd(i)=output_pass(1);
    deltaP_tube(i)=output_pass(2);
    deltaP_ehd(i)=output_pass(3);
    P_in_tube(i)=output_pass(4);
    P_out_tube(i)=output_pass(5);
    P_ehd(i)=P_out_tube(i)-deltaP_ehd(i);
    Q_t(i)=output_pass(6);
    T_ehd(i)=output_pass(7);
    h_ehd(i)=output_pass(8);
    T_out_air(i,:)=output_array(1,:);
    deltaP_a(i,:)=output_array(2,:);
    error(i)=2.*((deltaP_ihd(i)+deltaP_tube(i))-(deltaP_ehd(i-1)+deltaP_tube(i-1)))./ ...
        ((deltaP_ihd(i)+deltaP_tube(i))+(deltaP_ehd(i-1)+deltaP_tube(i-1)));   
%     disp('G_tube0=')
%     disp(G_tube0')
end
error(1)=(sum(G_tube0)-G_pass)/G_pass;
disp('error=')
disp(error)
% output_tube(1,:)=G_tube0;
% output_tube(2,:)=Q_t;
% output_tube(3,:)=deltaP_ihd;
% output_tube(4,:)=deltaP_ehd;
% output_tube(5,:)=P_in_tube;
% output_tube(6,:)=P_ihd;
% output_tube(7,:)=P_ehd;
% output_tube(8,:)=T_ehd;
% output_tube(9,:)=deltaP_tube;
end
