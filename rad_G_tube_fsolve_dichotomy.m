function [output_tube, output_result_distibution]=rad_G_tube_fsolve_dichotomy...
    (G_tube0,G_pass, tube_pass, P_inpass0, T_inlet0, va, ta, z, geometry)
%% basic information
% 该函数是基于压力平衡得到各个管路的流量分布
% 采用二分法迭代计算真实流量值
% 初值的选取很重要，最好是预先计算几个简单的看下，这里用的是fsolve函数求解后的结果，至于fsolve的结果是否可信，未知
% 需要提醒的是fsolve函数的求解速度是和变量的数量有关系的，变量越多，方程越复杂，求解难度也越大。
% 调用rad_ihd_tube_ehd单管程计算子函数
% 输出数据组计算微元分布（output_result_distibution）顺序说明：
% 1-换热量、2-热水压降、3-出口压力（热水）、4-空气对流换热系数、5-空气出口温度
% 6-重力压降、7-Nu数、8-出口温度（热水）、9-换热量误差、10-导热系数
% 11-传热系数、12-NTU数、13-最小热容、14-epsilon数、15-空气压降、16-摩擦压降、17-f因子
% Coded by Bei.W
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
output_result_distibution=ones(tube_pass,z,17);
%% main calculation codes
% G_tube0=G_pass./(tube_pass.*ones(1,tube_pass))'; % 第一根管初始赋值 kg/s
error=ones(1,tube_pass); % 误差初始赋值
G_1_max=0.04; % 之前计算结果
G_1_min=0.01; %
iter_1=0;
while abs(error(1))>0.01
    iter_1=iter_1+1;
    if iter_1>20
        break
    end
    G_tube0(1)=(G_1_max+G_1_min)/2;
    mr_ihd(1)=G_pass;
    mr_ehd(1)=G_tube0(1);
    P_ihd(1)=P_inpass0;
    T_ehd0=0; % Temperature of oulet in first tube
    [output_pass,output_array] = rad_ihd_tube_ehd(mr_ihd(1), G_tube0(1), mr_ehd(1), ...
        P_ihd(1), T_inlet0, T_ehd0, va(1,:), ta(1,:), z, geometry);
    deltaP_ihd(1)=output_pass(1);
    deltaP_tube(1)=output_pass(2);
    deltaP_ehd(1)=output_pass(3);
    P_in_tube(1)=output_pass(4);
    P_out_tube(1)=output_pass(5);
    Q_t(1)=output_pass(6);
    T_ehd(1)=output_pass(7);
    P_ehd(1)=output_pass(8);
    for o=1:length(output_array(1,:))
        output_result_distibution(1,:,o)=output_array(:,o);
    end
    for i=2:tube_pass
        G_i_max=0.02;
        G_i_min=0.005;
        error(i)=1;
        iter_i=0;
        while abs(error(i))>0.001
            iter_i=iter_i+1;
            if iter_i>20
                break
            end
            G_tube0(i)=(G_i_max+G_i_min)/2;
            mr_ihd(i)=mr_ihd(i-1)-G_tube0(i);
            mr_ehd(i)=mr_ehd(i-1)+G_tube0(i);
            P_ihd(i)=P_in_tube(i-1)-deltaP_ihd(i-1);
            [output_pass,output_array] = rad_ihd_tube_ehd(mr_ihd(i), G_tube0(i), mr_ehd(i),...
                P_ihd(i), T_inlet0, T_ehd(i-1), va(i,:), ta(i,:), z, geometry);
            deltaP_ihd(i)=output_pass(1);
            deltaP_tube(i)=output_pass(2);
            deltaP_ehd(i)=output_pass(3);
            P_in_tube(i)=output_pass(4);
            P_out_tube(i)=output_pass(5);
            Q_t(i)=output_pass(6);
            T_ehd(i)=output_pass(7);
            P_ehd(i)=output_pass(8);
            for o=1:length(output_array(1,:))
                output_result_distibution(i,:,o)=output_array(:,o); % store data into 3-dimension array
            end
            error(i)=2.*((deltaP_ihd(i)+deltaP_tube(i))-(deltaP_ehd(i-1)+deltaP_tube(i-1)))./ ...
                ((deltaP_ihd(i)+deltaP_tube(i))+(deltaP_ehd(i-1)+deltaP_tube(i-1)));
            if error(i)>0
                G_i_max=G_tube0(i);
            else
                G_i_min=G_tube0(i);
            end
        end
        running_iter=['大循环第' num2str(iter_1) '次-小循环第' num2str(i) '次-迭代次数为' num2str(iter_i)];
        disp(running_iter)
    end
    error(1)=(sum(G_tube0)-G_pass)/G_pass;
    if error(1)>0
        G_1_max=G_tube0(1);
    else
        G_1_min=G_tube0(1);
    end
    disp('error=')
    disp(error')
end
%% output parameters
output_tube(1,:)=G_tube0;
output_tube(2,:)=Q_t;
output_tube(3,:)=deltaP_ihd;
output_tube(4,:)=deltaP_ehd;
output_tube(5,:)=P_in_tube;
output_tube(6,:)=P_ihd;
output_tube(7,:)=P_ehd;
output_tube(8,:)=T_ehd;
output_tube(9,:)=deltaP_tube;
end
