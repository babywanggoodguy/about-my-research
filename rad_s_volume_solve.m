function err_Q=rad_s_volume_solve(IniV,IT,IP,G_t,va,ta,geometry)
% 此函数为求解计算微元内的流动换热
% 输入参数：假定出口压力及微元换热量、进口状态压力、焓值、风速、风温
% 输出参数：换热量误差、压降误差
%顺序情况下的微元计算
% volume=[OT,OP]
% ../  ..*是指对应元素的运算
% 温度单位均为℃
%% "input heat-exchanger variables"散热器基本输入参数
%% 换热器参数
L=geometry(1);%606.*0.001;%单流程长度m
H=geometry(2);%换热器高度
z=geometry(3);%50;%单流程换热微元数
L_t=geometry(4);%L./z;%单个换热单元的长度m
Lp=geometry(5);%1.3.*0.001;   %"louver pitch"百叶窗间距m
Ll=geometry(6);%6.8.*0.001;   %"louver length"百叶窗长度m
deg=geometry(7);%27;  % "louver angle"百叶窗角度deg
Fp=geometry(8);%1.3.*0.001;  %"fin pitch"肋片间距m
Fd=geometry(9);%16.*0.001;   %"fin width"肋片宽度m
Fh=geometry(10);%8.*0.001;  %"fin height"肋片高度m
delta_f=geometry(11);%0.1.*0.001;  %"fin thickness"肋片厚度m
Td=geometry(12);%16.*0.001;  % "flattube width"单通道平行管宽度m
Th=geometry(13);%1.8.*0.001;  % "flattube height"单通道平行管高度m
Cd=geometry(14);%1.3.*0.001;  % "innerchannel width"微通道宽度m
Ch=geometry(15);%1.6.*0.001;  % "innerchannel height"微通道高度m
delta_t=geometry(16);%0.1.*0.001;% "tube thickness"平行管厚度m
n=geometry(18);%10;  % "innerchannel number./flattube"内通道数量m
tube_sum=geometry(19);%84;

%% "calculate geometry value of tube and fin"
A_tr=n.*Cd.*Ch;% 内通道面积
W=2.*(Cd+Ch).*n;
Ar=2.*n.*(Cd+Ch).*L_t;%制冷剂侧换热面积
Aba=2.*(Td+Th).*L_t;%空气侧平行管换热面积
Afa=2.*Fd.*Fh./Fp.*L_t;%空气侧翅片换热面积
Aa=Aba+Afa;%空气侧总换热面积
Lh=0.5.*Lp.*tan(deg./180.*3.141592654);
Dhr=4.*A_tr./W;%制冷剂侧水力直径
Dha=2.*(Fp-delta_f).*(Fh-delta_f)./((Fp-delta_f)+(Fh-delta_f));%空气侧水力直径
% alpha=Ch./delta_t;
T_a_in=ta;%空气进口温度[C]
Ba=geometry(20);%露点温度[C]
Cp_a=CoolProp.HAPropsSI('cp_ha','T',T_a_in+273.15,'P',101325,'R',Ra); %空气进口比热容
rou_a=1./CoolProp.HAPropsSI('Vha','T',T_a_in+273.15,'P',101325,'R',Ra); %空气进口密度
err_Q=1;
OT_0=IniV(1);
OT_a=IniV(2);
%% 计算主体
    T_unit=(OT_0+IT)./2;
    prop_T=Property_Ethylene_Glycol(T_unit); %物性参数
    rho_r=prop_T(1);
    Cp=prop_T(2);
    lambda_t=prop_T(3);
    vis=prop_T(4);
    Pr=vis.*rho_r.*Cp/lambda_t.*1000;
    v_max=va.*(Fp+delta_f).*(Fh+Th)./((Fp-Lh-delta_f).*(Fh-delta_f));
    %% 风量计算
    A=L.*H./z./(tube_sum);
    Ga_v=va.*A;
    G_a=Ga_v.*rou_a;
    %% "calculate the h_a and pressure drop on the air-side" 我觉得这里有点问题,是否应该假设空气出口温度
    OT_a=abs(OT_0-IT).*Cp.*G_t./(Cp_a.*G_a)+T_a_in; %微元空气出口温度，℃
    T_a_unit=(OT_a+T_a_in)./2;
    Prop_air=Property_Air(T_a_unit);
    rou_a_m=Prop_air(1);
    v_a_m=Prop_air(4);
    c_a_m=Prop_air(2);
    lamb_a_m=Prop_air(3);
    lambda_Al=237;    
    %% 风侧换热系数及压降计算
    Pr_a=c_a_m.*v_a_m./lamb_a_m;
    Re_a=v_max.*Dha.*rou_a_m./v_a_m;
    % j_a=0.249.*Re_a.^(-0.42).*Lh.^(0.33).*(Ll./Fh).^1.1.*Fh.^0.26; % Davenport 1983, 100<Re<4000
    % j_a=Re_a^(-0.49).*(deg./90)^0.27.*(Fp./Lp)^(-0.14).*(Fh./Lp)^(-0.29).*(Td./Lp)^(-0.23).*(Ll./Lp)^0.68.*((Fh+Th)./Lp)^(-0.28).*(delta_f./Lp)^(-0.05); %Wang&chang 1997
    j_a=0.26712.*Re_a^(-0.194).*(deg./90)^0.257.*(Fp./Lp)^(-0.5177).*(Fh./Lp)^(-1.9045).*(Td./Lp)^(-0.2147).*(Ll./Lp)^1.7159.*(delta_f./Lp)^(-0.05);%董军启博士论文
    h_a=1.*rou_a_m.*v_max.*c_a_m.*j_a.*Pr_a^(-2./3);    
    %% "air-side pressure drop
    %     sigma=((Fp-Lh-delta_f)*(Fh-delta_f))/((Fp+delta_f)*(Fh+Th));
    %     Ke=0.0393*sigma^3-0.4963*sigma^2+0.0637*sigma+0.395; % 突扩压力损失系数
    %     Kc=0.1001*sigma^3+0.7806*sigma^2-1.8651*sigma+0.9823; % 突缩压力损失系数
    if Re_a<=130
        fa_1=14.39*Re_a^((-0.805)*Fp/Fh)*(log(1+(Fp/Lp)))^3.04;
        fa_2=((log(0.9+(delta_f/Fp)^0.48))^(-1.435))*(Dha/Lp)^(-3.01)*(log(0.5*Re_a))^(-3.01);
        fa_3=(Fp/Ll)^(-0.308)*(Fd/Ll)^(-0.308)*exp(-0.1167*(Th+Fh)/Th)*deg^0.35;
        f_a=fa_1*fa_2*fa_3;
    elseif Re_a>=230
        fa_1=4.97*Re_a^(0.6049-1.064/deg^0.2)*(log(0.9+(delta_f/Fp)^0.5))^(-0.527);
        fa_2=((Dha/Lp)*log(0.3*Re_a))^(-2.966)*(Fp/Ll)^((-0.7931)*(Th*Fh)/Fh);
        fa_3=((Th+Fh)/Fh)^(-0.0446)*log(1.2+(Lp/Fp)^1.4)^(-3.553)*deg^(-0.477);
        f_a=fa_1*fa_2*fa_3;
    else
        w=3.6-0.02*Re_a;
        fa_1=14.39*130^((-0.805)*Fp/Fh)*(log(1+(Fp/Lp)))^3.04;
        fa_2=((log(0.9+(delta_f/Fp)^0.48))^(-1.435))*(Dha/Lp)^(-3.01)*(log(0.5*130))^(-3.01);
        fa_3=(Fp/Ll)^(-0.308)*(Fd/Ll)^(-0.308)*exp(-0.1167*(Th+Fh)/Th)*deg^0.35;
        f_a_130=fa_1*fa_2*fa_3;
        fa_1=4.97*230^(0.6049-1.064/deg^0.2)*(log(0.9+(delta_f/Fp)^0.5))^(-0.527);
        fa_2=((Dha/Lp)*log(0.3*230))^(-2.966)*(Fp/Ll)^((-0.7931)*(Th*Fh)/Fh);
        fa_3=((Th+Fh)/Fh)^(-0.0446)*log(1.2+(Lp/Fp)^1.4)^(-3.553)*deg^(-0.477);
        f_a_230=fa_1*fa_2*fa_3;
        f_a=sqrt(0.5*((1+w)*f_a_130^2+(1-w)*f_a_230^2));
    end
    %     deltaP_air=0.5*rou_a_1*v_max^2*(((Fp+delta_f)*(Fh+Th)*rou_a_1*f_a)/((Fp-Lh-delta_f)*(Fh-delta_f)*rou_a_m)+(Kc+1-sigma^2)+(1-sigma^2-Ke)); %(Pa) from Yu-Juei Chang and Chi-Chuan Wang 2006
    %     f_a=0.54486.*Re_a^(-0.3068).*(deg./90)^0.444.*(Fp./Lp)^(-0.9925).*(Fh./Lp)^0.5458.*(Td./Lp)^0.0688.*(Ll./Lp)^(-0.2003); %董军启博士论文
    deltaP_air=4.*f_a.*Td.*rou_a_m.*v_max^2./Dha;
    Re_ref=(abs(G_t)./A_tr).*Dhr./(vis);
    if Re_ref<=2300
        Nu_t=4.36;
    else
        ff=(1.82.*log10(Re_ref)-1.64)^(-2)./8;% 赵宇单相换热公式
        Nu_t=ff.*(Re_ref-1000).*Pr./(1+12.7.*ff^0.5.*(Pr^(2./3)-1)); % Gnielinski公式
        %Nu_t=ff.*0.5.*(Re-1000).*Pr./(1+12.7.*(ff^0.5.*(Pr^(2./3)-1))^(0.5)).*(1+(Dhr./L_t)^(2./3));
    end
    %Nu_t=0.85.*ff.*(Re-1000).*Pr./(1+12.7.*(ff.*(Pr^(2./3)-1))^(0.5)).*(1+(Dhr./L_t)^(2./3));
    %Nu_t=0.023.*Re^0.8.*Pr^0.3334;
    k=Nu_t.*lambda_t./Dhr;
    %"pressure drop in the single regime"
    %% 单相区摩擦系数公式1 (祁照岗、赵宇文章中使用）
    %         if (Re_ref<=2300)
    %             f1=64./Re_ref;
    %         else
    %             f1=0.3164.*Re_ref^(-0.25);
    %         end
    %% 单相区摩擦系数公式2（Blasius关联式）
    if (Re_ref<=2500)
        f1=16./Re_ref;
    elseif (Re_ref<=20000)&&(Re_ref>2500)
        f1=0.079.*Re_ref^(-0.25);
    else
        f1=0.046.*Re_ref^(-0.2);
    end
    DeltaP_fr=f1.*4.*(abs(G_t)./A_tr)^2.*L_t./(2.*rho_r.*Dhr).*0.001;
    %     DeltaP_G=L_t.*g.*rou_r.*0.001;% 竖直管考虑重力压降
    DeltaP_G=0; % 水平管不考虑重力压降
    DeltaP_ref=DeltaP_fr+DeltaP_G; % 制冷剂侧压降 kPa
    OP_cal=IP-DeltaP_ref; %微元冷却液出口压力
    %% Churchill单相压降公式
    % C=(7./Re)^0.9+0.27.*0.000005./Dhr;
    % A=(2.457.*log(1./C))^16;
    % B=(37530./Re)^16;
    % f=2.*((8./Re)^12+(A+B)^(-1.5))^(1./12);
    % deltaP=f.* L_t.*(abs(G_t)./A_tr)^2./(2.*rou_r.*Dhr).*0.001;
    %% e-NTU 方法
    U=1./(Aa./(k.*Ar)+delta_t.*Aa./(lambda_Al.*Ar)+1./h_a);
    y1=abs(G_t).*Cp;
    y2=G_a.*c_a_m.*0.001;
    y_min=min(y1,y2);
    y_max=max(y1,y2);
    C_r=y_min./y_max;
    NTU=U.*Aa.*0.001./y_min;
    epsilon=1-exp(NTU^0.22.*(exp(-C_r.*NTU^0.78)-1)./C_r); % for single phase
    Qc=abs(epsilon.*y_min.*(IT-T_a_in));%e-NTU方法计算所得换热量 /Kw
    Q0=abs(G_t.*Cp.*(OT_0-IT)); % 假设微元换热量
    err_Q=2.*(Q0-Qc)./(Q0+Qc) % 换热量计算误差
    if err_Q > 0
        OT_min=OT_0;
    else
        OT_max=OT_0;
    end
end
%% output parameters
volume_result(1)=Q0;
volume_result(2)=DeltaP_ref;
volume_result(3)=OP_cal;
volume_result(4)=h_a;
volume_result(5)=OT_a;
volume_result(6)=DeltaP_G;
volume_result(7)=Nu_t;
volume_result(8)=OT_mid;
volume_result(10)=k;
volume_result(11)=U;
volume_result(12)=NTU;
volume_result(13)=y_min;
volume_result(14)=epsilon;
volume_result(15)=deltaP_air;
volume_result(16)=DeltaP_fr;
volume_result(17)=f_a;
end