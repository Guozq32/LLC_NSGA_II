clear all
clc
close all
% format long
global Po Vo Vin_min Vin_max rou miu ds diso Rdson alpha beta kcof etac etam pop_size pm V M xl xu Ae Aw Ve AeAw h hw Icm Sm DL tdead Coss eta tf Ron Izvs Qrr
% Icm利兹线电流；Sm利兹线股数；DL 利兹线外径；

%%  将稳态工作点的计算数据设置为全局变量
global Vin_given  fs_min fs_max Cr Ln Z0 Z1 wr wn ILs_rms_max Ipeak_max ...
ILs_rms_Matrix_less Ipeak_t0t2_Matrix_less Ipeak_t2t3_Matrix_less Ipeak_Matrix_less fs_Matrix_less  ...
r0_Matrix_less phi0_Matrix_less the0_Matrix_less  ...
r1_Matrix_less phi1_Matrix_less the1_Matrix_less M_Matrix_less Vin_Matrix_less  ...
ILs_rms_Matrix_more Ipeak_Matrix_more fs_Matrix_more  r0_Matrix_more  phi0_Matrix_more the0_Matrix_more  ...
r1_Matrix_more phi1_Matrix_more the1_Matrix_more M_Matrix_more Vin_Matrix_more

global M_max M_min

etac = 0.30;                  % distribution index for crossover 交叉因子
etam = 0.3;                  % distribution index for mutation / mutation constant %% 突变的分布指数
M=2; %优化目标个数
pop_size=200;           % Population size 种群数为偶数
gen_max=50;            % MAx number of generations - stopping criteria 产生种群的代，即执行遗传算法次数
Vo=200;
Vin_min=150;
Vin_max=250;
% n=Vo./Vm;
Po=1000;
rou=2.12e-8;	%%铜的电阻率单位0.0185Ω·mm2/m->@25℃;  0.0212Ω·mm2/m->@70℃;
miu=4*pi*1e-7; %%真空磁导率
tdead=60e-9;
Qrr= 20e-9;   %%% 二极管反向恢复电荷，通过手册查找
Coss=260e-12;
Izvs=2*Vin_max*Coss/tdead;%满足软开关电流限制值
eta=1;


% 关断损耗
Qgd_p=7e-9;                        %栅漏电荷                   单位C  
Qgs_p=8e-9;                        %栅源电荷                   单位C  
Ciss_p=640e-12;                    %输入电容                   单位F  

Crss_p=2.3e-12;                      %逆导电容                   单位F  
Vsp_p=6;                            %米勒平台电压（估计值）     单位V  
VGS_th_p=2.8;                       %栅极阈值电压               单位V  
RG_int_p=6;                       %内部栅极电阻               单位Ω

Vcc=15;                             %驱动电压正电压        单位V  
Vee=-5;                             %驱动电压负电压        单位V  
Rdriver_off=2;                    %驱动电阻（估计值）    单位Ω
RG_off_p=Rdriver_off+RG_int_p;      %关断电阻
Cgs_p=Ciss_p-Crss_p;                %内部栅源电容
Qth_off_p=(VGS_th_p-Vee)*Cgs_p;
QG_off_p=Qgd_p+Qgs_p-Qth_off_p;
IG_off_p=(Vsp_p-Vee)/RG_off_p;
tf=QG_off_p/IG_off_p;           %关断时电流和电压的交越时间



ds=0.1e-3;  %%利兹线每根直径 单位m
diso=0.1e-3; %%每层绕组之间绝缘距离
Ron=0.08;
alpha=1.658963082164934; %%这应该是TDKPC47的参数，国产的PC系列看损耗图其他参数一样大概会是1.5-2倍的损耗。
ki=0.123934484240003;
beta=2.410283969308211;
fun=@(theta) power(abs(cos(theta)),alpha)*power(2,beta-alpha);
q=integral(fun,0,2*pi);
kcof=ki/(power(2*pi,alpha-1)*q);


%     1      2       3       4       5       6
%  PQ26/25  PQ32/20 PQ32/30 PQ35/35 PQ40/40 PQ50/50
Ae = [120	169     167     190     202     328]*1e-6;    %%单位m^2   磁芯截面积
Aw = [84.53	80.50	149.1	220     325.98	433.20]*1e-6; %%单位m^2   窗口面积
AeAw=[10143 13604.5 24899.7 41800   65846.95 142089.6]*1e-12; %%单位m^2   AP
Ve = [6530  9440    12500   16300   20500   37100]*1e-9;  %%单位m^3 磁芯体积
h  = 0.7*[16.1  11.5    21.3    25      29.5    36.1]*1e-3; %单位m 窗口高度 有效高度利用率取70%
hw = 0.7*[5.3   7       7       8.8     11.1    12]*1e-3;   %单位m 窗口宽度 有效宽度利用率取70%
%     1     2   3    4   5    6   7     8    9    10    11  12  13
Icm= [3.95 4.3  4.7  5.1 5.5  5.9 6.3   7.05 7.85 8.65  9.8	11	11.8]; %利兹线对应电流有效值 分别是100到300股利兹线的耐流值
Sm = [100  110  120  130 140  150 160   180  200  220   250	280 300]; %每根利兹线股数
DL = [1.4  1.47 1.53 1.6 1.66 1.71 1.77 1.88 1.98 2.08  2.21 2.34 2.42]*1e-3; %单位m 对应线的外径

%   1   2      3     4      5     6       %   1 2   3  4    5    6    
xl=[0.7 8e-6  60e3  2  0.13 0.13];   %%x=[n,Lr,fr,K,  Bm_T,Bm_L]  
xu=[1.5 50e-6 140e3  7.5  0.25 0.25];    %%x=[n,Lr,fr,K,  Bm_T,Bm_L]
V=size(xl,2)
pm=1/V;                     % Mutation Probability 变异率



%%                  1                   2                         3                      4                5                 6                  7                     8               9                10                  11         12          13             14              15              16         17   18   19      20         21        22      23        24             25                 26                           27          28          29        30            31         32          33          34    35      36       37
%% Para=     [      N1                  N2                        N3                     Nl1              Nl2              Nl3                 m1                    m2              m3              Rac1               Rac2        Rac3        Sm1             Sm2             Sm3            Km_T        Km_L Ve_T Ve_L    Tr_index Lr_index    Lm      Cr     Ploss_avg  Pswitch_off_avg P_conduction_Mosfet_avg P_conduction_Diode_avg P_Tr_Cu_avg P_Lr_Cu_avg P_Tr_Fe_avg P_Lr_Fe_avg  Ipeak_max ILs_rms_max  fs_max fs_min  Range_fs error];
%% Para_cell=[  Ploss_Cell       Pswitch_off_Cell     P_conduction_Mosfet_Cell P_conduction_Diode_Cell  P_Tr_Cu_Cell  P_Lr_Cu_Cell    P_Tr_Fe_t0t2_Cell      P_Tr_Fe_t2t3_Cell   P_Lr_Fe_t0t2_Cell P_Lr_Fe_t2t3_Cell ILs_rms_Cell  Im1_Cell    fs_Cell       Vin_Cell     Vin_Matrix1_Cell Bm_T_t0t2_Cell];
xl_temp=repmat(xl, pop_size,1);  %% repmat([1 2; 3 4],2,3) 返回一个矩阵。xl_temp是pop_size×V维数的矩阵  行数代表种群数，列数代表优化目标变量数
xu_temp=repmat(xu, pop_size,1);
x = xl_temp+((xu_temp-xl_temp).*rand(pop_size,V)); %%rand(pop_size,V)是pop_size×V维数的矩阵；x仍然是pop_size×V维数的矩阵
[x] = Re_random(x,xl,xu) %%就是确保x的每个量的，判断准谐振结束时重新分配励磁电感两端电压值大小，当接近跨入PN模态时重新随机当前组参数

%% Evaluate objective function
[x,Para,Para_cell,ff,err] = Cal_Loss(x);%ff 两列优化目标，err违反值


error_norm=normalisation(err);                     % Normalisation of the constraint violation  %%生成的error_norm矩阵为pop_size×1 维的矩阵，把err矩阵所有元素除以每一列的最大值，再把每一行所有元素相加
population_init=[x ff error_norm];
%[population front]=NDS_CD_cons(population_init);    % Non domination Sorting on initial population 对初始种群进行非支配排序
[population,initial_data_Para_sorted,initial_data_Para_cell_sorted,front]=NDS_CD_cons(population_init,Para,Para_cell);    % Non domination Sorting on initial population 对初始种群进行非支配排序
population_initial=population;
population_front1=population(find(population(:,V+M+2)==1),:);
population_1_feasible=population(find(population(:,V+M+1)==0),:);
population_1_infeasible=population(find(population(:,V+M+1)~=0),:);
parent_data_Para=initial_data_Para_sorted;
parent_data_Para_cell=initial_data_Para_cell_sorted;
%% Generation Starts
for gen_count=1:gen_max
    gen_count
% selection (Parent Pt of 'N' pop size)   %% population矩阵是排完序的样本，其中还包括误差、层级、和拥挤度
parent_selected=tour_selection(population);                     % 10 Tournament selection 锦标赛竞争法  %%随机挑选出交叉变异的父样本， 父样本数与population维数相同
% Reproduction (Offspring Qt of 'N' pop size)
child_offspring  = genetic_operator(parent_selected(:,1:V));    % SBX crossover and polynomial mutation 模拟二进制交叉和多项式变异； %%输入的向量只包含优化的变量，不包含优化目标函数、误差、层级和拥挤度；
%子代维数与随机生成的父代维数相同
[xc,Parac,Parac_cell,ffc,errc] = Cal_Loss(child_offspring);
child_data_Parac=Parac;
child_data_Parac_cell=Parac_cell;
error_norm=normalisation(errc);  %%子代误差标幺化                                
child_offspring=[xc ffc error_norm]; %% 生成包含目标值和标幺化误差的子代矩阵

%% INtermediate population (Rt= Pt U Qt of 2N size)
data_Para_inter=[parent_data_Para; child_data_Parac]; %%拼接父代子代的数据
data_Para_cell_inter=[parent_data_Para_cell ; child_data_Parac_cell]; %%拼接父代子代的数据
population_inter=[population(:,1:V+M+1) ; child_offspring(:,1:V+M+1)]; %% V+1列到V+M列对应对应ff和fff的列数(优化结果目标)，V+M+1对应err的列数
% [population_inter_sorted,front]=NDS_CD_cons(population_inter);              % Non domination Sorting on offspring  对父代和子代组成的矩阵重新进行非支配排序
[population_inter_sorted,data_Para_sorted,data_Para_cell_sorted,front]=NDS_CD_cons(population_inter,data_Para_inter,data_Para_cell_inter);  
%% Replacement - N 根据分层筛选下一代昂本
% new_pop=replacement(population_inter_sorted, front);
[new_pop,data_Para_new,data_Para_cell_new]=replacement(population_inter_sorted,data_Para_sorted,data_Para_cell_sorted, front);
population=real(new_pop);
parent_data_Para=data_Para_new;
parent_data_Para_cell=data_Para_cell_new;

 if(gen_count==2)
    population_front3=population(find(population(:,V+M+2)==1),:);
    population_3_feasible=population(find(population(:,V+M+1)==0),:);
    population_3_infeasible=population(find(population(:,V+M+1)~=0),:);
    figure (2)
    plot(population_front3(:,V+1),population_front3(:,V+2),'*');
    elseif(gen_count==4)  %第5代  
    population_front5=population(find(population(:,V+M+2)==1),:);
    population_5_feasible=population(find(population(:,V+M+1)==0),:);
    population_5_infeasible=population(find(population(:,V+M+1)~=0),:);
    elseif (gen_count==9)  %第10代
    population_front10=population(find(population(:,V+M+2)==1),:);
    population_10_feasible=population(find(population(:,V+M+1)==0),:);
    population_10_infeasible=population(find(population(:,V+M+1)~=0),:);
    elseif (gen_count==14)  %第15代
    population_front15=population(find(population(:,V+M+2)==1),:);
    population_15_feasible=population(find(population(:,V+M+1)==0),:);
    population_15_infeasible=population(find(population(:,V+M+1)~=0),:);  
    elseif (gen_count==19)  %第20代
    population_front20=population(find(population(:,V+M+2)==1),:);
    population_20_feasible=population(find(population(:,V+M+1)==0),:);
    population_20_infeasible=population(find(population(:,V+M+1)~=0),:);            
    elseif (gen_count==24)  %第25代
    population_front25=population(find(population(:,V+M+2)==1),:);
    population_25_feasible=population(find(population(:,V+M+1)==0),:);
    population_25_infeasible=population(find(population(:,V+M+1)~=0),:);  
    elseif (gen_count==29)  %第30代
    population_front30=population(find(population(:,V+M+2)==1),:);
    population_30_feasible=population(find(population(:,V+M+1)==0),:);
    population_30_infeasible=population(find(population(:,V+M+1)~=0),:);            
    elseif (gen_count==34)  %第35代
    population_front35=population(find(population(:,V+M+2)==1),:);
    population_35_feasible=population(find(population(:,V+M+1)==0),:);
    population_35_infeasible=population(find(population(:,V+M+1)~=0),:);
    elseif (gen_count==39)  %第40代
    population_front40=population(find(population(:,V+M+2)==1),:);
    population_40_feasible=population(find(population(:,V+M+1)==0),:);
    population_40_infeasible=population(find(population(:,V+M+1)~=0),:);
    elseif (gen_count==44)  %第45代
    population_front45=population(find(population(:,V+M+2)==1),:);
    population_45_feasible=population(find(population(:,V+M+1)==0),:);
    population_45_infeasible=population(find(population(:,V+M+1)~=0),:);   
   elseif (gen_count==49)  %第50代
    population_front50=population(find(population(:,V+M+2)==1),:);
    population_50_feasible=population(find(population(:,V+M+1)==0),:);
    population_50_infeasible=population(find(population(:,V+M+1)~=0),:);
end

end

[new_pop,sort]=sortrows(new_pop,V+1); %% sortrows(new_pop,V+1) 按照第V+1列的升序对矩阵new_pop进行重新排列； %% V+1对应第一个目标值，按照第一个目标值进行排序
Parar=parent_data_Para(sort,:); %% 按sort的顺序同样进行数据的排序
Parar_cell=parent_data_Para_cell(sort,:); %% 按sort的顺序同样进行数据的排序
% new_pop=sortrows(new_pop,V+1); %% sortrows(new_pop,V+1) 按照第V+1列的升序对矩阵new_pop进行重新排列； %% V+1对应第一个目标值，按照第一个目标值进行排序
% [xr,Parar,Parar_cell,ffr,errr] = Cal_Loss(population);
Result_pop=new_pop(:,1:9);
[Result_population,Result_data_Parar_sorted,Result_data_Parar_cell_sorted,Result_front]=NDS_CD_cons(Result_pop, Parar, Parar_cell); 
Result_pop_feasible=Result_pop(find(population(:,V+M+1)==0),:);
Result_pop_infeasible=Result_pop(find(population(:,V+M+1)~=0),:);
flag=find(Result_pop(:,V+M+1)==0);
% xr=xc;
% Parar=Parac;
% Parar_cell=Parac_cell;
% ffr=ffc;
% errr=errc;
% Result_pop=[xr ffr errr];
% [Result_population,Result_front]=NDS_CD_cons(Result_pop); 

for i=1:size(Result_front(1).fr,2)%Result_front(1).fr,2是pareto最优面的个数，rank=1；

Front_pop(i,:)=Result_pop(flag(Result_front(1).fr(i)),:);%Front_pop是 [pareto最优面的个数，V+M+1]
Front_Parac(i,:)=Parar(flag(Result_front(1).fr(i)),:); %Front_pop是 [pareto最优面的个数，Parac]
Front_Parac_cell(i,:)=Parar_cell(flag(Result_front(1).fr(i)),:); %Front_pop是 [pareto最优面的个数，Parac]
end
% for i=1:size(Result_front(1).fr,2)
% Front_pop(i,:)=Result_pop(Result_front(1).fr(i),:);
% Front_Parar(i,:)=Parar(Result_front(1).fr(i),:);
% Front_Parar_cell(i,:)=Parar_cell(Result_front(1).fr(i),:);
% end
% [Front_pop,Front_pop_sort]=sortrows(Front_pop,V+1);
% Front_Parar=Front_Parar(Front_pop_sort,:);
% Front_Parar_cell=Front_Parar_cell(Front_pop_sort,:);

% for i=1:size(front(1).fr,2)
% Front_pop_l(i,:)=population_inter(front(1).fr(i),:);
% end
%% Result and Pareto plot
figure(1),plot(Result_pop(:,V+1),Result_pop(:,V+2),'*'); hold on %最后一次所有迭代出的最后一代100个 
figure(1),plot(Front_pop(:,V+1),Front_pop(:,V+2),'ro'); hold on %pareto最优面
figure(2),plot(Front_pop(:,V+1),Front_pop(:,V+2),'*'); hold on %横坐标是损耗，纵坐标是体积

figure(2),plot(population_front40(:,V+1),population_front40(:,V+2),'*'); hold on %横坐标是损耗，纵坐标是体积

