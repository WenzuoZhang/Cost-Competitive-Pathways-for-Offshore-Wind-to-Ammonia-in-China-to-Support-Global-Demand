clear clc
yalmip('clear')

%管道布局优化
load inxy3.mat
load data_v_NH.mat
load TypicalWeek.mat
load TypicalWeekidx.mat
final_result=[];
y=inxy3(:,1);%纬度
x=inxy3(:,2);%经度
Num_max = 89 + (ceil(111.1*cosd(x(:))*0.25*0.5)-13)*5;

R=[1.000 	0.872 	0.796 	0.749 	0.722 	0.706 	0.697 	0.694 	0.693 
1.000 	0.815 	0.710 	0.649 	0.613 	0.593 	0.582 	0.578 	0.577 
1.000 	0.761 	0.633 	0.562 	0.521 	0.498 	0.486 	0.481 	0.481 
1.000 	0.710 	0.565 	0.486 	0.443 	0.418 	0.406 	0.401 	0.400 
1.000 	0.663 	0.504 	0.421 	0.376 	0.352 	0.339 	0.334 	0.333 
1.000 	0.620 	0.449 	0.364 	0.319 	0.295 	0.283 	0.278 	0.278];%BAU
% R=[1.000 	0.761 	0.633 	0.562 	0.521 	0.498 	0.486 	0.481 	0.481 
% 1.000 	0.710 	0.565 	0.486 	0.443 	0.418 	0.406 	0.401 	0.400 
% 1.000 	0.663 	0.504 	0.421 	0.376 	0.352 	0.339 	0.334 	0.333 
% 1.000 	0.620 	0.449 	0.364 	0.319 	0.295 	0.283 	0.278 	0.278 
% 1 	0.579 	0.401 	0.316 	0.271 	0.248 	0.237 	0.232 	0.231 
% 1 	0.540 	0.358 	0.273 	0.231 	0.208 	0.198 	0.193 	0.192];%OTA
t=0;
BAUorOTL = 1;%1BAU 2OTA
Ele_effi=[0.65 0.68 0.71 0.74 0.77 0.80;0.65 0.7 0.75 0.8 0.85 0.9];
Ele_lowbound=[0.2 0.18 0.16 0.14 0.12 0.1;0.2 0.16 0.12 0.08 0.04 0];
Ele_upbound=[1 1.03 1.06 1.09 1.12 1.15;1 1.06 1.12 1.18 1.24 1.30];
Ele_ramp=[0.5 0.6 0.7 0.8 0.9 1;0.5 0.65 0.8 1 1 1];
HB_effi=[0.82 0.83 0.84 0.85 0.86 0.87;0.82 0.84 0.86 0.88 0.9 0.92];
HB_lowbound=[0.3 0.27 0.24 0.21 0.18 0.15;0.3 0.25 0.2 0.15 0.1 0.05];
HB_upbound=[1 1.03 1.06 1.09 1.12 1.15;1 1.05 1.1 1.15 1.2 1.35];
HB_ramp=[0.3 0.35 0.4 0.45 0.5 0.55;0.3 0.4 0.5 0.6 0.7 0.8];
for wangge = 1:2565
   wangge
   tic
P_wind_ty=P_wind_input(:,:,wangge);
WT_Rated_power = Num_max(wangge)*12/1000; %风电装机容量 10^6kW
idx=idx_input(:,wangge);
Rated_speed=11; %当前风机型号额定风速
Cut_in_speed=3; %切入风速
Cut_out_speed=25; %切出风速
Frictioncoeffic=1/7; %地表摩擦系数
Heightofturbine=200; %海上风机轮毂高度200米
Heightoftower=100;  %海上测风塔高度100米，因此需要将风速递推至200米
ty=inxy3(wangge,1);%当前风电场的经度
tx=inxy3(wangge,2);%当前风电场的纬度
ti=89-(tx-3)/0.25;%对应位置
tj=1+(ty-100)/0.25;%对应位置
Wind_speed = squeeze(data_v_NH(ti,tj,:));%100米测风塔风速
Net_Wind_speed = Wind_speed.*(Heightofturbine/Heightoftower)^(Frictioncoeffic); %风速高度转换为200米
for i=1:length(Wind_speed) %风机实际功率计算 10^6kW
    if Net_Wind_speed(i)<Cut_in_speed || Net_Wind_speed(i)>Cut_out_speed
        P_wind(i,1)=0;
    elseif Net_Wind_speed(i)>=Cut_in_speed && Net_Wind_speed(i)<=Rated_speed
        P_wind(i,1)=WT_Rated_power.*((Net_Wind_speed(i)-Cut_in_speed)./(Rated_speed-Cut_in_speed));
    elseif Net_Wind_speed(i)>Rated_speed && Net_Wind_speed(i)<=Cut_out_speed
        P_wind(i,1)=WT_Rated_power;
    end
end

D = pip_layout_prim(Num_max(wangge)); %最角上风机和网格中心的距离
Cost_intrayarray=D*5230000/10^6*R(3,t/5+1);%10^6元

%聚类可视化
for i = 1:52
    idx_day((i-1)*7+1:i*7,1) = repmat(idx(i),7,1);
    day_idx((i-1)*7+1:i*7,1) = (idx(i)-1)*7+[1:7];
end
idx_day(365,1)=idx_day(358,1);
day_idx(365,1)=day_idx(358,1);

idx_month(1,1:31) = idx_day(1:31,1);
idx_month(2,1:28) = idx_day(32:59,1);
idx_month(3,1:31) = idx_day(60:90,1);
idx_month(4,1:30) = idx_day(91:120,1);
idx_month(5,1:31) = idx_day(121:151,1);
idx_month(6,1:30) = idx_day(152:181,1);
idx_month(7,1:31) = idx_day(182:212,1);
idx_month(8,1:31) = idx_day(213:243,1);
idx_month(9,1:30) = idx_day(244:273,1);
idx_month(10,1:31) = idx_day(274:304,1);
idx_month(11,1:30) = idx_day(305:334,1);
idx_month(12,1:31) = idx_day(335:365,1);

%风机成本计算
Cost_Plant=(1300+1184)*6.76*R(3,t/5+1);
Total_Cost_Plant=Cost_Plant*WT_Rated_power; %当前网格所需风机成本（10^6元）
Resvalue_rate=0.05; %贴现率
Net_Cost_Plant=Total_Cost_Plant*(1-Resvalue_rate); %当前风电场购置风机净支出成本
Self_capa_rate=0.2; %业主投资比例
Life=15; % 折旧15年
Depreciation_charge=[ones(1,15)*Net_Cost_Plant/Life,zeros(1,5)]; %年折旧（直线折旧）
Maintenance_cost=[zeros(1,5) ,ones(1,5)*Total_Cost_Plant*0.012,ones(1,10)*Total_Cost_Plant*0.015]; %年维修成本
Wages_and_benefits=ones(1,20)*36*60000*(1+0.14+0.31+0.12)/10^6; %员工年薪  10^6元
Material=ones(1,20)*WT_Rated_power*20; %材料费
Insurance=ones(1,20)*Total_Cost_Plant*0.0025; %保险（保票）费
Other_cost=ones(1,20)*WT_Rated_power*30; %其他成本
Load_rate=0.05; %偿债利率
Load_z=Total_Cost_Plant*(1-Self_capa_rate); %负债总量
Annual_Changhuaibenjin=Load_z/15; %年偿还本金
Construction_Load=Load_z*ones(1,15)-Annual_Changhuaibenjin*[0:1:14]; %运营期年对应负债量
Interest=Construction_Load*Load_rate; %年利息
Benchmark_rate=0.08; %折现率，基准利润率
Total_Cost=Depreciation_charge+Maintenance_cost+Wages_and_benefits+Material+Insurance+Other_cost+[Interest zeros(1,5)]; %年总成本
Operating_Cost=Total_Cost-Depreciation_charge-[Interest zeros(1,5)];
Resvalue=[zeros(1,15),-Total_Cost_Plant*Resvalue_rate,zeros(1,4)];
Total_Cost_Plant_Self=Total_Cost_Plant*Self_capa_rate;
Jiekuanbenjinchanghuan_AA_z=[ones(1,15)*Annual_Changhuaibenjin,zeros(1,5)];%年借款本金偿还
Initial_investment=Total_Cost_Plant_Self+sum(Jiekuanbenjinchanghuan_AA_z./((1+Benchmark_rate).^[1:20])) ;
O_M=sum(Operating_Cost./((1+Benchmark_rate).^[1:20]));
Interest_total=sum([Interest zeros(1,5)]./((1+Benchmark_rate).^[1:20]));
Resvalue_total=sum(Resvalue./((1+Benchmark_rate).^[1:20])) ;
Total_Plant=Initial_investment+O_M+Interest_total+Resvalue_total ; %风电场总成本

%运行优化
time_long = length(P_wind_ty); 
yita_BT = 0.97; %蓄电池 充放效率
gama_BT = 0.001; %蓄电池 每小时损失0.1%
yita_HE = 0.95; %储氢罐储能效率为95%
gama_HE = 0.001/24; %储氢罐每天损失0.1%;

z = sdpvar(1,1,'full'); %目标函数分母的替换，实际意义为氨厂全生命周期产氨量之和的倒数
P_wind_sum=0;
for j = 1:4*7
    P_wind_sum = P_wind_sum + sum(P_wind_ty(24*(j-1)+1:24*j)*sum(day_idx==j,'all'));
end
z_true = P_wind_sum/( 1/(1000/176.5)/(Ele_effi(BAUorOTL,t/5+1)*0.0252) + 237.3/HB_effi(BAUorOTL,t/5+1)/1000)*9.8181; % z为 氨厂全生命周期之和(单位转换为10^10kg=1董^_^)的倒数

M_z =100; %一个恒大于z的极大值

Inst_Capa_Ele_z = sdpvar(1,1,'full'); %电解槽容量*z
P_Ele_z = sdpvar(time_long,1,'full'); %电解槽耗电功率*z 10^6kW
H_Ele_z = P_Ele_z*Ele_effi(BAUorOTL,t/5+1)*0.0252;  %电解槽产氢速率*z, 10^6kg
C = [Inst_Capa_Ele_z>=0]; %电解槽容量*z大于等于0 10^6kW
C = [C,Inst_Capa_Ele_z<=WT_Rated_power*z]; %电解槽容量*z小于风电装机总容量*z

index_Ele = binvar(time_long,1,'full'); %电解槽耗运行状态变量
index_Inst_z_Ele = sdpvar(time_long,1,'full'); %中间变量 电解槽耗运行状态变量*装机容量*z
M_Inst_Capa_Ele_z = 100*WT_Rated_power; %一个恒大于电解槽容量*z的极大值

% 电解槽耗运行状态变量 * （装机容量*z）的线性化过程
C = [C,index_Inst_z_Ele>=0];
C = [C,index_Inst_z_Ele<=M_Inst_Capa_Ele_z*index_Ele];
C = [C,index_Inst_z_Ele<=Inst_Capa_Ele_z];
C = [C,index_Inst_z_Ele>=Inst_Capa_Ele_z-M_Inst_Capa_Ele_z*(1-index_Ele)];
% 
C = [C,P_Ele_z>=Ele_lowbound(BAUorOTL,t/5+1)*index_Inst_z_Ele];%电解槽出力*z范围为0-1.6
C = [C,P_Ele_z<=Ele_upbound(BAUorOTL,t/5+1)*index_Inst_z_Ele];

% 制氢厂爬坡约束
for i = 1:4
    time_idx = (i-1)*7*24;
    for j = 1:7*24
        if j ==1
            C = [C, P_Ele_z(time_idx+j)-P_Ele_z(time_idx+7*24)<=Inst_Capa_Ele_z*Ele_ramp(BAUorOTL,t/5+1)];
            C = [C, P_Ele_z(time_idx+j)-P_Ele_z(time_idx+7*24)>=-Inst_Capa_Ele_z*Ele_ramp(BAUorOTL,t/5+1)];
        else
            C = [C, P_Ele_z(time_idx+j)-P_Ele_z(time_idx+j-1)<=Inst_Capa_Ele_z*Ele_ramp(BAUorOTL,t/5+1)];
            C = [C, P_Ele_z(time_idx+j)-P_Ele_z(time_idx+j-1)>=-Inst_Capa_Ele_z*Ele_ramp(BAUorOTL,t/5+1)];  
        end
    end
end

Inst_Capa_BT_z = sdpvar(1,1,'full'); %电池容量*z
P_BT_C_z = sdpvar(time_long,1,'full'); %电池储电功率*z
P_BT_D_z = sdpvar(time_long,1,'full'); %电池储电功率*z
P_BT_z = P_BT_D_z - P_BT_C_z; %电池功率*z

S_BT_z = sdpvar(time_long,1,'full'); %电池逐时储能量*z

C = [C,Inst_Capa_BT_z>=0]; %电池容量*z大于等于0 10^6kWh
C = [C,P_BT_C_z>=0];
C = [C,P_BT_D_z>=0];
C = [C,P_BT_C_z<=Inst_Capa_BT_z*0.5]; %电池供电功率*z范围为-0.1容量*z~0.1容量*z
C = [C,P_BT_D_z<=Inst_Capa_BT_z*0.5]; %电池供电功率*z范围为-0.1容量*z~0.1容量*z

for i = 1:4
    for j =1:7*24
        time_day = (i-1)*7*24 + j;
        if j == 1
            C = [C,S_BT_z(time_day,1)==0.5*Inst_Capa_BT_z*(1-gama_BT)-P_BT_D_z(time_day,1)/yita_BT+P_BT_C_z(time_day,1)*yita_BT] ;
        else
            C = [C,S_BT_z(time_day,1)==S_BT_z(time_day-1,1)*(1-gama_BT)-P_BT_D_z(time_day,1)/yita_BT+P_BT_C_z(time_day,1)*yita_BT] ;
        end
    end
end

C = [C,S_BT_z>=0.1*Inst_Capa_BT_z];%储能量*z范围为0.1装机容量*z-0.9装机容量*z
C = [C,S_BT_z<=0.9*Inst_Capa_BT_z];

C = [C,S_BT_z(1*7*24)>=0.5*Inst_Capa_BT_z];%储能量末大于等于初
C = [C,S_BT_z(2*7*24)>=0.5*Inst_Capa_BT_z];%储能量末大于等于初
C = [C,S_BT_z(3*7*24)>=0.5*Inst_Capa_BT_z];%储能量末大于等于初
C = [C,S_BT_z(4*7*24)>=0.5*Inst_Capa_BT_z];%储能量末大于等于初

%%
Inst_Capa_HE_z = sdpvar(1,1,'full'); %储氢罐容量*z 10^6kg
H_HE_C_z = sdpvar(time_long,1,'full'); %储氢罐储氢速率*z 10^6kg
H_HE_D_z = sdpvar(time_long,1,'full'); %储氢罐供氢速率*z 10^6kg
H_HE_z = H_HE_D_z - H_HE_C_z; %储氢罐逐时储氢量*z 10^6kg
S_HE_z =sdpvar(time_long,1,'full'); 

C = [C,Inst_Capa_HE_z>=0]; %储氢罐容量*z大于等于0
C = [C,H_HE_C_z>=0];
C = [C,H_HE_D_z>=0];
C = [C,H_HE_C_z<=Inst_Capa_HE_z*0.5]; %储氢罐容量-功率系数为0.1
C = [C,H_HE_D_z<=Inst_Capa_HE_z*0.5]; %储氢罐供氢速率*z范围为-0.1容量*z~0.1容量*z

for i = 1:4
    for j =1:7*24
        time_day = (i-1)*7*24 + j;
        if j == 1
            C = [C,S_HE_z(time_day,1)==0.5*Inst_Capa_HE_z*(1-gama_HE)-H_HE_D_z(time_day,1)/yita_HE+H_HE_C_z(time_day,1)*yita_HE] ;
        else
            C = [C,S_HE_z(time_day,1)==S_HE_z(time_day-1,1)*(1-gama_HE)-H_HE_D_z(time_day,1)/yita_HE+H_HE_C_z(time_day,1)*yita_HE] ;
        end
    end
end

C = [C,S_HE_z>=0.1*Inst_Capa_HE_z]; %储能量范围*z为0.1装机容量*z-0.9装机容量*z
C = [C,S_HE_z<=0.9*Inst_Capa_HE_z];

C = [C,S_HE_z(1*7*24)>=0.5*Inst_Capa_HE_z];%储能量末大于等于初
C = [C,S_HE_z(2*7*24)>=0.5*Inst_Capa_HE_z];%储能量末大于等于初
C = [C,S_HE_z(3*7*24)>=0.5*Inst_Capa_HE_z];%储能量末大于等于初
C = [C,S_HE_z(4*7*24)>=0.5*Inst_Capa_HE_z];%储能量末大于等于初

Inst_Capa_AA_z = sdpvar(1,1,'full'); %制氨厂装机容量*z 
A_AA_z = sdpvar(time_long,1,'full'); %制氨厂产氨速率*z 10^6kg/h
H_AA_z = A_AA_z/(1000/176.5); %制氨厂耗氢速率*z 10^6kg
P_AA_z = A_AA_z*237.3/HB_effi(BAUorOTL,t/5+1)/1000; %制氨厂耗电功率*z 10^6kW [http://tjj.nanjing.gov.cn/material/njnj_2019/gongye/7-13.htm]
C = [C,Inst_Capa_AA_z>=0]; %制氨厂装机容量*z大于等于0

index_AA = binvar(time_long,1,'full'); % Haber-Bosch与空气分离设备 运行状态变量
index_Inst_z_AA = sdpvar(time_long,1,'full'); %中间变量 制氨厂运行状态变量*装机容量*z
M_Inst_Capa_AA_z = 100*max(P_wind_ty)/(237.3/HB_effi(BAUorOTL,t/5+1)/1000); %一个恒大于制氨厂容量*z的极大值

% 制氨厂运行状态变量 * （装机容量*z）的线性化过程
C = [C,index_Inst_z_AA>=0];
C = [C,index_Inst_z_AA<=M_Inst_Capa_AA_z*index_AA];
C = [C,index_Inst_z_AA<=Inst_Capa_AA_z];
C = [C,index_Inst_z_AA>=Inst_Capa_AA_z-M_Inst_Capa_AA_z*(1-index_AA)];

C = [C,P_AA_z>=HB_lowbound(BAUorOTL,t/5+1)*index_Inst_z_AA]; %制氨速率范围为0.3-1装机容量
C = [C,P_AA_z<=HB_upbound(BAUorOTL,t/5+1)*index_Inst_z_AA];

for i = 2:time_long
    C = [C, P_AA_z(i)-P_AA_z(i-1)<=Inst_Capa_AA_z*HB_ramp(BAUorOTL,t/5+1)];
    C = [C, P_AA_z(i)-P_AA_z(i-1)>=-Inst_Capa_AA_z*HB_ramp(BAUorOTL,t/5+1)];
end

C = [C,P_Ele_z+P_AA_z<=P_BT_z+P_wind_ty*z]; %电 -- 能源供需平衡
C = [C,H_AA_z<=H_Ele_z+H_HE_z]; %氢 -- 能源供需平衡

A_AA_z_sum=0; %年累积产氢量
A_AA_z_sum_season_max = sdpvar(1,1,'full'); %各季节累积产氢量最大值
for j = 1:4*7
    C = [C,A_AA_z_sum_season_max>=sum(A_AA_z(24*(j-1)+1:24*j)*sum(day_idx==j,'all'))];
    A_AA_z_sum = A_AA_z_sum + sum(A_AA_z(24*(j-1)+1:24*j)*sum(day_idx==j,'all'));
end

C = [C,A_AA_z_sum*9.8181==z_true]; %z与A_AA的关系
C = [C,z>=0];
C = [C,z<=M_z];

%成本计算
% PEM电解槽成本计算
Cost_Ele=2000*R(1,t/5+1); %电解槽投资成本 800美元/10^6kW
Net_Cost_Ele_z=Inst_Capa_Ele_z*Cost_Ele*(1-Resvalue_rate);
Lifecycle=15; %寿命长度
Cost_water_year=0.005*sum(H_Ele_z)*9*52/4/10^6*R(1,t/5+1);%水处理费用
Cost_hycomp=7348*7.7*sum(H_Ele_z)*52/4/10^6*R(4,t/5+1);%氢压缩费用
%工商细价
Depreciation_charge2_z=[ones(1,15)*Net_Cost_Ele_z/Lifecycle zeros(1,5)];
Maintenance_cost2_z=[zeros(1,5) ones(1,5)*Inst_Capa_Ele_z*Cost_Ele*0.012 ones(1,5)*Inst_Capa_Ele_z*Cost_Ele*0.015 zeros(1,5)];
Wages_and_benefits2=ones(1,20)*12*60000/10^6 *(1+0.14+0.31+0.12); %常数 千元
Material2_z=ones(1,20)*Inst_Capa_Ele_z*20;
Insurance2_z=ones(1,20)*Net_Cost_Ele_z*0.0025;
Other_cost2_z=ones(1,20)*Inst_Capa_Ele_z*30;
Load_z=Inst_Capa_Ele_z*Cost_Ele*(1-Self_capa_rate);%负债
Annual_Changhuaibenjin2_z=Load_z/15;
Construction_Load2_z=Load_z*ones(1,15)-Annual_Changhuaibenjin2_z*[0:1:14];
Interest2_z=Construction_Load2_z*Load_rate;
Total_Cost2_z=Depreciation_charge2_z+Maintenance_cost2_z+Wages_and_benefits2*z+Material2_z+Insurance2_z+Other_cost2_z+[Interest2_z zeros(1,5)]; %鎬绘垚鏈垂鐢?
Operating_Cost2_z=Total_Cost2_z-Depreciation_charge2_z-[Interest2_z zeros(1,5)];
Resvalue2_z=[zeros(1,15),-Inst_Capa_Ele_z*Cost_Ele*Resvalue_rate zeros(1,4)];

Balance_construction_z=Inst_Capa_Ele_z*Cost_Ele*Self_capa_rate;
Repayment_selfowned_z=Balance_construction_z/Lifecycle;

s_z=[];
for i=1:20
    S_z=Balance_construction_z-Repayment_selfowned_z;
    Balance_construction_z=S_z;
    s_z=[s_z S_z];
end
End_balance_z=s_z;
Beginning_balance_z=[Inst_Capa_Ele_z*Cost_Ele*Self_capa_rate s_z];
Benefit_self_z=Beginning_balance_z(1:20)*0.08;

Income_tax2_z=[zeros(1,3) Benefit_self_z(4:6)/(1-0.25/2)-Benefit_self_z(4:6) Benefit_self_z(7:20)/(1-0.25)-Benefit_self_z(7:20)];
VAT_input_tax2_z=(Maintenance_cost2_z+Material2_z+Insurance2_z+Other_cost2_z)/(1+0.07)*0.07;
VAT_output_tax2_z=(Total_Cost2_z+Benefit_self_z)/(1+0.07)*0.07;
VAT_tax2_z=VAT_output_tax2_z-VAT_input_tax2_z;
Urban_tax2_z=VAT_tax2_z*0.05;
Education_tax2_z=VAT_tax2_z*0.05;

Total_Cost_Ele_Self_z=Inst_Capa_Ele_z*Cost_Ele*Self_capa_rate;
Jiekuanbenjinchanghuan2_z=[ones(1,15)*Annual_Changhuaibenjin2_z,zeros(1,5)];
Tax2_z=Income_tax2_z+VAT_tax2_z+Urban_tax2_z+Education_tax2_z;

Initial_investment2_z=Total_Cost_Ele_Self_z+sum(Jiekuanbenjinchanghuan2_z./((1+Benchmark_rate).^[1:20])) ;
O_M2_z=sum(Operating_Cost2_z./((1+Benchmark_rate).^[1:20])) ;
Interest_total2_z=sum([Interest2_z zeros(1,5)]./((1+Benchmark_rate).^[1:20])) ;
Resvalue_total2_z=sum(Resvalue2_z./((1+Benchmark_rate).^[1:20])) ;
Tax_Total2_z=sum(Tax2_z./((1+Benchmark_rate).^[1:20]));
Replace_cost_z=Inst_Capa_Ele_z*Cost_Ele*0.8/(1+Benchmark_rate)^15;
Total_Ele_z=Initial_investment2_z+O_M2_z+Interest_total2_z+Resvalue_total2_z+Tax_Total2_z+Replace_cost_z;

% 氨厂成本计算
Cost_AA=4709.090909 * 6.76*R(1,t/5+1); %AA投资成本 4709.090909美元/(10^6kg/h) [A feasibility study of implementing an Ammonia Economy]
Total_AA_Plant_z=Inst_Capa_AA_z*Cost_AA; %当前网格所需AA成本（千元）
Net_Cost_AA_z=Total_AA_Plant_z*(1-Resvalue_rate); %购置AA净支出成本
Life_AA=20; %折旧时长
%工商细价
Depreciation_charge_AA_z=ones(1,20)*Net_Cost_AA_z/Life_AA; %年折旧成本（直线折旧）
Maintenance_cost_AA_z=[zeros(1,5) ones(1,5)*Total_AA_Plant_z*0.012 ones(1,5)*Total_AA_Plant_z*0.015 zeros(1,5)];
Wages_and_benefits_AA=ones(1,20)*12*60000/10^6*(1+0.14+0.31+0.12); %员工年薪
Material_AA_z=ones(1,20)*Inst_Capa_AA_z*20; %材料费
Insurance_AA_z=ones(1,20)*Net_Cost_AA_z*0.0025; %保险（保票）费
Other_cost_AA_z=ones(1,20)*Inst_Capa_AA_z*30; %其他成本
Load_AA_z=Total_AA_Plant_z*(1-Self_capa_rate); %负债总量
Annual_Changhuanbenjin_AA_z=Load_AA_z/Life_AA; %折旧期内年偿还本金
Construction_Load_AA_z=Load_AA_z*ones(1,Life_AA)-Annual_Changhuanbenjin_AA_z*[0:1:Life_AA-1]; %运营期各年剩余负债量
Interest_AA_z=Construction_Load_AA_z*Load_rate; %各年所需偿还的负债利息
Total_Cost_AA_z=Depreciation_charge_AA_z+Maintenance_cost_AA_z+Wages_and_benefits_AA*z+Material_AA_z+Insurance_AA_z+Other_cost_AA_z+Interest_AA_z; %年成本
Operating_Cost_AA_z=Total_Cost_AA_z-Depreciation_charge_AA_z-Interest_AA_z; % 年实际运营成本
Resvalue_AA_z=[zeros(1,15),-Total_AA_Plant_z*Resvalue_rate zeros(1,4)]; %补贴偿还

Total_Cost_AA_Self_z=Total_AA_Plant_z*Self_capa_rate; %业主实际投资成本
Jiekuanbenjinchanghuan_AA_z=ones(1,Life_AA)*Annual_Changhuanbenjin_AA_z;%年借款本金偿还
Initial_investment_AA_z=Total_Cost_AA_Self_z+sum(Jiekuanbenjinchanghuan_AA_z./((1+Benchmark_rate).^[1:20])) ;%业主总投资量贴现值
O_M_AA_z=sum(Operating_Cost_AA_z./((1+Benchmark_rate).^[1:20])); %运行成本贴现值
Interest_total_AA_z=sum(Interest_AA_z./((1+Benchmark_rate).^[1:20])); %负债利息贴现值
Resvalue_total_AA_z=sum(Resvalue_AA_z./((1+Benchmark_rate).^[1:20])) ; %补贴偿还贴现值
Total_AA_z=Initial_investment_AA_z+O_M_AA_z+Interest_total_AA_z+Resvalue_total_AA_z ; %AA总投资贴现值

%电池成本计算
Total_BT_z = Inst_Capa_BT_z * 200*R(1,t/5+1) + Inst_Capa_BT_z * 0.1 * 1500*R(1,t/5+1); %200/kwh,1500元/kw xueqiu.com/9804362435/194616922
%储氢罐成本计算
Total_HE_z = Inst_Capa_HE_z * 7000*R(1,t/5+1); %7000*(10^6元/10^6kg) [Ammonia to power: Forecasting the levelized cost of electricity from green ammonia in large-scale power plants]
%船舶运输成本计算
Total_B_z = 0.35; %0.35元/Kg氨气
%海上氨厂平台成本
Total_AF = 6.76*10000*10000/10^6*1.5*R(2,t/5+1); %
%低温储氨容器成本
Total_AS_z = A_AA_z_sum_season_max * 20.2*10^6*6.76/25000/1000*R(3,t/5+1);
%场内电缆成本
Total_intrayarray = Cost_intrayarray;

%水处理
Total_water = sum(Cost_water_year./((1+Benchmark_rate).^[1:20]));
%压缩
Total_hycomp = sum(Cost_hycomp*ones(1,20)./((1+Benchmark_rate).^[1:20]));

p_curt = sum(- P_Ele_z - P_AA_z+P_BT_z+P_wind_ty*z);
H_curt = sum(-H_AA_z+H_Ele_z+H_HE_z);

%优化目标
obj = (Total_BT_z + Total_HE_z + Total_AA_z + Total_Plant*z +Total_Ele_z + Total_intrayarray*z +Total_AF*z +Total_AS_z +Total_water+Total_hycomp)/z_true; %制氨LOCA
ops = sdpsettings('solver' ,'gurobi');
ops.savesolveroutput =1;
ops.SolutionLimit=1;
ops.gurobi.MIPGap = 0.05;
result1 = optimize(C,obj,ops);%调用cplex参数设置

if result1.problem ~=12
    f_min = value(obj);%得到最优解，则输出OBJ的值
     mip_gap = result1.solveroutput.result.mipgap;
else
    f_min = -Inf;%未得到最优解，则输出Inf，无穷大
    mip_gap = Inf;
end
toc
z =value(z);
Inst_Capa_Ele = value(Inst_Capa_Ele_z)/z; %电解槽装机容量 Kw
Inst_Capa_BT = value(Inst_Capa_BT_z)/z;
Inst_Capa_HE = value(Inst_Capa_HE_z)/z;
Inst_Capa_AA = value(Inst_Capa_AA_z)/z;
Total_power = sum(P_wind,'all'); %年总发电量 Kwh
Total_hour_use = Total_power/WT_Rated_power; %年风电等效利用小时数 h
Total_A = z_true/z/1000; %全寿命周期总制氨量 t
Total_A_year = Total_A/9.8181; %年总制氨量 t

H_Ele = value(H_Ele_z/z); %电解槽典型周出力
H_EL_all = zeros(8760,1);
opt_day_index = linspace(1,length(day_idx),length(day_idx));
for i = 1:length(opt_day_index)
    ty_day = day_idx(i,1); %当前自然日对应的典型日编号
    ty_week = floor(ty_day/7)+1;%当前自然日对应的典型周编号
    starttime_ty_day = 7*24*(ty_week-1)+24*(mod(ty_day,7)-1)+1; %典型周编号对应的起止出力序列
    if i == 1
        H_EL_all(1,1)=H_Ele(starttime_ty_day,1) ;
        for j =2:24
            H_EL_all(j,1)=H_EL_all(j-1,1)+H_Ele(starttime_ty_day+j-1,1) ;
        end
    else
        starttime_day = (i-1)*24+1; %当前自然日的起止时间 如第2个自然日为 25-48，第3个自然日为 73-96小时
        for j = 1:24
            H_EL_all(starttime_day+j-1,1) = H_EL_all(starttime_day+j-2,1) + H_Ele(starttime_ty_day+j-1,1) ;
        end
    end
end
Total_H_year = H_EL_all(end)/1000;%年总制氢量  t
Total_H = Total_H_year*9.8181;%全寿命周期总制氢量 t

Total_Ele = value(Total_Ele_z)/z/10000; %电解槽建设成本 万元
Total_AA = value(Total_AA_z)/z/10000; %氨厂建设成本 万元
Total_BT = value(Total_BT_z)/z/10000; %电储能建设成本 万元
Total_HE = value(Total_HE_z)/z/10000; %储氢罐建设成本 万元
Total_B = value(Total_B_z)/z/10000; %船舶总费用 万元
Total_AS = value(Total_AS_z)/z/10000; %储氨总费用 万元

Loca_wind = Total_Plant*z/z_true; %LOCA 风电场分项
Loca_Ele = value(Total_Ele_z)/z_true; %LOCA 电解槽分项
Loca_AF = Total_AF*z/z_true; %LOCA 氨厂平台分项
Loca_AA = value(Total_AA_z)/z_true; %LOCA 氨厂建设分项
Loca_pip = Total_intrayarray*z/z_true; %LOCA 海底管道分项
Loca_BT = value(Total_BT_z)/z_true; %LOCA 电储能分项
Loca_HE = value(Total_HE_z)/z_true; %LOCA 储氢罐分项
Loca_B = Total_B_z/z_true; %LOCA 船舶分项
Loca_AS = value(Total_AS_z)/z_true;%LOCA 储氨分项
Lcoa_water = value(Total_water)/z_true;%水处理分项
Lcoa_hycomp = value(Total_hycomp)/z_true;%氢压缩分项

LCOH=(Loca_wind+Loca_Ele+Loca_AF+Loca_pip+Loca_BT+Loca_HE+Lcoa_water+Lcoa_hycomp)*5.66;

LOCA = f_min; %LOCA 制1kg氨成本
toc
Result = [wangge,ty,tx,WT_Rated_power,Inst_Capa_Ele,Inst_Capa_BT,Inst_Capa_HE,Inst_Capa_AA,Total_power,Total_hour_use,Total_H,Total_H_year,Total_A,Total_A_year,...
    Total_Plant/10000,Total_Ele,Total_AF/10000,Total_AA,Total_intrayarray/10000,Total_BT,Total_HE,Total_B,Total_AS,...
    Loca_wind,  Loca_Ele, Loca_AF, Loca_AA, Loca_pip, Loca_BT, Loca_HE, Loca_B,Loca_AS,Lcoa_water,Lcoa_hycomp,...
    LOCA,LCOH,mip_gap,z];
final_result = [final_result;Result];
end