function [output_tube, output_result_distibution]=rad_G_tube_calc(G_tube0,G_pass, tube_pass, P_inpass0, T_inlet0, va, ta, z, geometry)
%% basic information
% �ú����ǻ���ѹ��ƽ��õ�������·�������ֲ�
% ���ö��ַ�����������ʵ����ֵ
% ����rad_ihd_tube_ehd���̼ܳ����Ӻ���
% ������������΢Ԫ�ֲ���output_result_distibution��˳��˵����
% 1-��������2-��ˮѹ����3-����ѹ������ˮ����4-������������ϵ����5-���������¶ȡ�6-����ѹ����7-Nu����8-�����¶ȣ���ˮ����9-��������10-����ϵ��
% 11-����ϵ����12-NTU����13-��С���ݡ�14-epsilon����15-����ѹ����16-Ħ��ѹ����17-f����
%% Ԥ�ȷ����ڴ棬����߼����ٶ�
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
error=ones(1,tube_pass); % ����ʼ��ֵ
output_result_distibution=ones(tube_pass,z,17);
%% main calculation codes
% G_tube0=G_pass./(tube_pass.*ones(1,tube_pass))'; % ��һ���ܳ�ʼ��ֵ kg/s
mr_ihd(1)=G_pass;
mr_ehd(1)=G_tube0(1);
P_ihd(1)=P_inpass0;
T_ehd0=0; % Temperature of oulet in first tube
[output_pass,output_array] = rad_ihd_tube_ehd(mr_ihd(1), G_tube0(1), mr_ehd(1), P_ihd(1), T_inlet0, T_ehd0, va(1,:), ta(1,:), z, geometry);
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
    mr_ihd(i)=mr_ihd(i-1)-G_tube0(i);
    mr_ehd(i)=mr_ehd(i-1)+G_tube0(i);
    P_ihd(i)=P_in_tube(i-1)-deltaP_ihd(i-1);
    [output_pass,output_array] = rad_ihd_tube_ehd(mr_ihd(i), G_tube0(i), mr_ehd(i), P_ihd(i), T_inlet0, T_ehd(i-1), va(i,:), ta(i,:), z, geometry);
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
    %         running_iter=['��ѭ��-' num2str(iter_1) '��-Сѭ��-' num2str(i) '��-��������-' num2str(iter_i)];
    %         disp(running_iter)
end
error(1)=(sum(G_tube0)-G_pass)/G_pass;
disp('error=')
disp(error')
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
