

function [Flag_error,error_1,Flag_PN,Flag_OPO,Flag_PON,Flag_NOP,Flag_ZVS,Flag_M_less,Flag_M_more ]= NSGA_LLC_steady_cal_fs_all_range(n,Lr,fr,Lm,R)


global pop_size Vo Vin_max Vin_min Izvs Coss tdead Po
global Vin_given fs_min fs_max Cr Ln Z0 Z1 wr wn ILs_rms_max Ipeak_max ...
ILs_rms_Matrix_less Ipeak_t0t2_Matrix_less Ipeak_t2t3_Matrix_less Ipeak_Matrix_less fs_Matrix_less  ...
r0_Matrix_less phi0_Matrix_less the0_Matrix_less  ...
r1_Matrix_less phi1_Matrix_less the1_Matrix_less M_Matrix_less Vin_Matrix_less  ...
ILs_rms_Matrix_more Ipeak_Matrix_more fs_Matrix_more  r0_Matrix_more  phi0_Matrix_more the0_Matrix_more  ...
r1_Matrix_more phi1_Matrix_more the1_Matrix_more M_Matrix_more Vin_Matrix_more

global M_max  M_min 


 M_max = Vo./(n*Vin_min);

 M_min = Vo./(n*Vin_max);

Vin_given=[150 160 170 180 190 200 210 220 230 240 250];


%%%%%%%%%%%%%(pop_size,1)不随开关频率改变%%%%%%%%%%%%%


Ln=zeros(pop_size,1);
Cr=zeros(pop_size,1);
Z0=zeros(pop_size,1);
Z1=zeros(pop_size,1);
wr0=zeros(pop_size,1);
wr1=zeros(pop_size,1);
wr=zeros(pop_size,1);
wn=zeros(pop_size,1);
In=zeros(pop_size,1);
Irm_estimate=zeros(pop_size,1);
Ir0=zeros(pop_size,1);
Vcr0=zeros(pop_size,1);
Ir0N=zeros(pop_size,1);
Vcr0N=zeros(pop_size,1);
Ir2=zeros(pop_size,1);
Vcr2=zeros(pop_size,1);
Ir2N=zeros(pop_size,1);
Vcr2N=zeros(pop_size,1);
t2=zeros(pop_size,1);
tr1=zeros(pop_size,1);
M=zeros(pop_size,1);


r0=zeros(pop_size,1);
the0=zeros(pop_size,1);
r1=zeros(pop_size,1);
the1=zeros(pop_size,1);
phi0=zeros(pop_size,1);
phi1=zeros(pop_size,1);
% VLm_fs_1=zeros(pop_size,1);

Vin=zeros(pop_size,1);

iteration_times=zeros(pop_size,1);
Range_min=zeros(pop_size,1);
R_opo=zeros(pop_size,1);
Vlm_pn=zeros(pop_size,1);
error_1=zeros(pop_size,1);
ILs_rms_max=zeros(pop_size,1);
Ipeak_max=zeros(pop_size,1);
fs_min=zeros(pop_size,1);
fs_max=zeros(pop_size,1);
x=zeros(7,pop_size);  %%% 7行 pop_size列数

%%%%%%  判断标志位的定义    %%%%%%
Flag_newton=zeros(pop_size,1);  %%% 牛顿迭代法出错标志位

Flag_M_less=zeros(pop_size,1);  %%% 满足最大增益M（欠谐振）要求的标志位
Flag_M_more=zeros(pop_size,1);  %%% 满足最小增益M（过谐振）要求的标志位
Flag_M=zeros(pop_size,1);       %%% 满足增益M要求的标志位，满足为1，初始化为0   



Flag_PN=zeros(pop_size,1);
Flag_OPO=zeros(pop_size,1);
Flag_PON=zeros(pop_size,1);     %%% 判断会不会从O模态进入N模态，即会不会进入PON模态，进入置1
Flag_NOP=zeros(pop_size,1);

Flag_error=zeros(pop_size,1);   %%% 表示该种群参数是否是满足条件的可行解，1不满足，0满足
Flag_ZVS =zeros(pop_size,1);        %%% 来表示

error_1 =zeros(pop_size,1);     %%%表示当下种群的不可行解违反值



m=zeros(pop_size,3);

r0_Matrix_less=zeros(pop_size,11);
phi0_Matrix_less=zeros(pop_size,11);
the0_Matrix_less=zeros(pop_size,11);
r1_Matrix_less=zeros(pop_size,11);
phi1_Matrix_less=zeros(pop_size,11);
the1_Matrix_less=zeros(pop_size,11);
M_Matrix_less=zeros(pop_size,11);
Vin_Matrix_less=zeros(pop_size,11);
fs_Matrix_less=zeros(pop_size,11);

r0_Matrix_more=zeros(pop_size,11);
phi0_Matrix_more=zeros(pop_size,11);
the0_Matrix_more=zeros(pop_size,11);
r1_Matrix_more=zeros(pop_size,11);
phi1_Matrix_more=zeros(pop_size,11);
the1_Matrix_more=zeros(pop_size,11);
M_Matrix_more=zeros(pop_size,11);
Vin_Matrix_more=zeros(pop_size,11);
fs_Matrix_more=zeros(pop_size,11);

%在扩展匝比n的范围后，会出现部分个体在输入电压范围内完全工作在PO或者NP模态，为了标记出这些个体，定义变量FLag_range
Flag_range=zeros(pop_size,2);% 如果第一位为1则表示全工作在PO模态（欠谐真），如果第二位为1，则表示全工作在NP模态（过谐振）


for i=1:pop_size
    
    %%在后续的损耗计算中，以10V为间隔取11个点计算，首先定义m，用于标记每个个体PO模态电压点的范围
    %%和NP模态电压点的范围，m(i,1-3)初始值相同，即PO和NP模态分界的电压点，第一列在欠谐振过程中被改变，第二列在过谐振过程中被改变，第三列不变
    if n(i)<=Vo/Vin_max  % 判断主要工作模态的范围，是否会出现全部PO模态或者NP模态
        m(i,1)=11;
        m(i,2)=11;
        m(i,3)=11;
        Flag_range(i,1)=1;%全工作在PO模态，所计算的11个点均为PO模态下
        Flag_M_more(i,1)=1;%如果全工作在欠谐振下，则最小增益要求不再考虑，只要满足最大的增益要求即为可行个体
    else if n(i)>=Vo/Vin_min
            m(i,1)=0;
            m(i,2)=0;
            m(i,3)=0;
            Flag_range(i,2)=1;%全工作在NP模态，所计算的11个点均为NP模态下
            Flag_M_less(i,1)=1;%如果全工作在过谐振下，则最大增益要求不再考虑，只要满足最小的增益要求即为可行个体
        else
   for j=1:1:10
    if ((n(i)*Vin_given(j))<=Vo)&&((n(i)*Vin_given(j+1))>Vo)
        m(i,1)=j;
        m(i,2)=j;
        m(i,3)=j;
    
    end
    end
        end
    end  
end
1
%% 通过循环对每个种群进行稳态值求解

for i=1:pop_size    

s1=0;
s2=0;


%%%%%%  根据传回的参数，计算第i种群下系统的电路参数值  %%%%%%

Vin(i)=Vo/n(i);
Cr(i)=1/4/pi/pi/fr(i)/fr(i)/Lr(i);
Ln(i)=Lm(i)/Lr(i);
Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i));wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
wr(i)=wr0(i);
wn(i)=wr1(i);
In(i)=Vin(i)/Z0(i);




%%%%%% 迭代前先排除OPO模态
R_OPO(i)=n(i)*n(i)*Vo*pi*Z0(i)/2/(Vo*(Lm(i)+Lr(i))/Lm(i)-n(i)*Vin(i)) ;

if R>=R_OPO(i)            %%% 可以留一定的裕量
    Flag_OPO(i)=1;
    error_1(i)=10;
    Flag_error(i)=1;
    continue                %%% 该种群参数不合适，退出该种群迭代
end

%%%%%% 迭代前PN模态

Vm=Lm(i)/(Lm(i)+Lr(i))*(2*R*Vin(i)-n(i)*Vo*pi*Z0(i))/2/R;
if Vm <= -Vo/n(i)
    Flag_PN(i)=1;
    error_1(i)=100;
    Flag_error(i)=1;
    continue                %%% 该种群参数不合适，退出该种群迭代
end

%% 牛顿迭代法

%%%%%%  在迭代求解稳态值前，根据传的参数估计准谐振状态下的迭代初始值  %%%%%%

Irm_estimate(i)=abs(Vo/Lm(i)/4/fr(i)/n(i));
Ir0(i)=-Irm_estimate(i);
Vcr0(i)=(Vin(i)-Vo/n(i)-n(i)*pi*Po*Z0(i)/2/Vo);%%vin-vo/n-npiPoZ0/2vo
Ir0N(i)=Ir0(i)/In(i);
Vcr0N(i)=Vcr0(i)/Vin(i);
Ir2(i)=Irm_estimate(i);
Vcr2(i)=-(Vin(i)-Vo/n(i)-n(i)*pi*Po*Z0(i)/2/Vo);
Ir2N(i)=Ir2(i)/In(i);
Vcr2N(i)=Vcr2(i)/Vin(i);
t2(i)=0.5/fr(i);
tr1(i)=0;

%%%%%牛顿迭代法参数设计%%%%%%%
eps1=1e-4;
eps2=1e-4;
lemda=0.5;

%% 欠谐振--牛顿迭代法
%%%%%%  把r0等待求解变量与各个初始值联系起来,设置迭代初始值     %%%%%%

%%%%%%  欠谐振，不同频率范围表达式不一样    %%%%%%   
M(i)=1; 
r0(i)=sqrt(Ir0N(i)^2+(Vcr0N(i)-(1-M(i)))^2);
the0(i)=atan(-Ir0N(i)/(Vcr0N(i)-(1-M(i))));
r1(i)=sqrt((1+Ln(i))*Ir2N(i)^2+(Vcr2N(i)-1)^2);
the1(i)=atan(-sqrt(1+Ln(i))*Ir2N(i)/(Vcr2N(i)-1));
phi0(i)=t2(i)*wr0(i);
phi1(i)=tr1(i)*wr1(i);   

%定义七个未知待求解参数r0,phi0,the0,r1,phi1,the1,M
x(1:7,i)=[r0(i),phi0(i),the0(i),r1(i),phi1(i),the1(i),M(i)]; % 7xPopsize

step=-500;
j=1;

for j_fs=(fr(i)):step:0
    %%判断匝比n，在输入电压范围内是否可存在PO模态，用m(i)判断
    if m(i,3)==0
        break %%表示需要记录的11个电压点均在NP模态下，跳出欠谐振牛顿迭代
    end
    
    %%%%%随着扫频改变周期%%%%%%%
    Ts(i,j)=1/j_fs;
    fs_Matrix1_less(i,j)=j_fs;
    iteration_times(i)=iteration_times(i)+1;
%     if j_fs==(fr(i)-500)
%         x(2,i)=pi+0.01;
%         x(5,i)=0.005;
%     end

        for k=1:5000
        
        b=F_less(x(:,i));%b输出是7个方程，方程理论上均=0

        norm_b=norm(b,inf);%矩阵的无穷范数是–各元素先取绝对值而后按行相加的最大值 ，本质上就是取最大的abs(A)，得到一个最大绝对值。
        if norm_b<eps1 %输出的7个方程中的最大值小于约定值则视为收敛完成跳出循环
            break;
        end
        A=Jac_less(x(:,i)); %A是导数
        dx=-pinv(A)*b; %见牛顿下山法公式(本文lamda没有下山逐渐减小，这样的优点是更快缺点是容易不收敛)
        x(:,i)=x(:,i)+lemda*dx; %新的迭代参数
        norm_dx=norm(dx,inf);%矩阵的无穷范数是–各元素先取绝对值而后按行相加的最大值
        if norm_dx<eps2
            break;
        end
        end


    %%%%%%  将迭代后的数据记录在矩阵当中      %%%%%%

    r0_Matrix1_less(i,j)=x(1,i);
    phi0_Matrix1_less(i,j)=x(2,i);
    the0_Matrix1_less(i,j)=x(3,i);
    r1_Matrix1_less(i,j)=x(4,i);
    phi1_Matrix1_less(i,j)=x(5,i);
    the1_Matrix1_less(i,j)=x(6,i);
    M_Matrix1_less(i,j)=x(7,i);    
    Vin_Matrix1_less(i,j)=Vo/n(i)/M_Matrix1_less(i,j);  %%%%%% 可能没必要

    if k==5000
            Flag_newton(i)=1; %Flag_1=1代表收敛失败
            error_1(i)=100;
            Flag_error(i)=1;
            sa(i)=1;
            break
    end

if abs(M_Matrix1_less(i,1)-1)>5e-2 || (r0_Matrix1_less(i,j)<0) || (r1_Matrix1_less(i,j)<0 )    %%%说明准谐振下增益不为1，迭代错误
    error_1(i)=100;
    Flag_error(i)=1;
    sb(i)=1;
    break
end


    %%%%%%  判断是否满足软开关条件     %%%%%%
            if abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*sin(the0_Matrix1_less(i,j))/Z0(i))<(2*Vin_Matrix1_less(i,j)*Coss/tdead)
%                 Flag_error(i)=1;
                Flag_ZVS(i)=1;
                error_1(i)=1;
                break
            end

    %%%%%%      判断是否进入PON模态，即当前迭代结果是否可靠     %%%%%%
    
    %%%%%计算励磁电感电压最小值%%%%%%%    
    if (the1_Matrix1_less(i,j)<pi && (phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))>pi)
        Vlm_min(i,j)=-Lm(i)*Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)/(Lm(i)+Lr(i));
    else
        Vlm_min(i,j)=min(Lm(i)*Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*cos(the1_Matrix1_less(i,j))/(Lm(i)+Lr(i)),Lm(i)*Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*cos(phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))/(Lm(i)+Lr(i)));  
    end 

    %%%%%%  根据O模态下的最小励磁电感电压来判断会不会从O模态进入N模态
    if Vlm_min(i,j)<=-Vo/n(i)

        Flag_PON(i)=1;
        error_1(i) =1;
        Flag_error(i)=1;
        break
    end

    %%%%%计算电流值%%%%%%%        ------包括输出电流平均值，谐振电流有效值以及谐振电流峰值
%     Io_avg_Matrix1_less(i,j)=2*Vin_Matrix1_less(i,j)*(r0_Matrix1_less(i,j)*cos(the0_Matrix1_less(i,j))-r0_Matrix1_less(i,j)*cos(phi0_Matrix1_less(i,j)+the0_Matrix1_less(i,j))-r0_Matrix1_less(i,j)*sin(the0_Matrix1_less(i,j))*phi0_Matrix1_less(i,j)-M_Matrix1_less(i,j)*Lr(i)*power(phi0_Matrix1_less(i,j),2)/2/Lm(i))/n(i)/wr(i)/Ts(i,j)/Z0(i);
    ILs_rms_Matrix1_less(i,j)=sqrt((r0_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*(2*phi0_Matrix1_less(i,j)+sin(2*the0_Matrix1_less(i,j))-sin(2*phi0_Matrix1_less(i,j)+2*the0_Matrix1_less(i,j)))/wr(i)/Ts(i,j)/Z0(i)/Z0(i)/2)+(r1_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*(2*phi1_Matrix1_less(i,j)+sin(2*the1_Matrix1_less(i,j))-sin(2*phi1_Matrix1_less(i,j)+2*the1_Matrix1_less(i,j)))/wn(i)/Ts(i,j)/Z1(i)/Z1(i)/2));
    if (the0_Matrix1_less(i,j)<pi/2 && (phi0_Matrix1_less(i,j)+the0_Matrix1_less(i,j))>pi/2)||((the0_Matrix1_less(i,j)<-pi/2 && (phi0_Matrix1_less(i,j)+the0_Matrix1_less(i,j))>-pi/2))
        Ipeak_t0t2_Matrix1_less(i,j)=abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)/Z0(i));
    else
        Ipeak_t0t2_Matrix1_less(i,j)=max(abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*sin(the0_Matrix1_less(i,j))/Z0(i)),abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*sin(phi0_Matrix1_less(i,j)+the0_Matrix1_less(i,j))/Z0(i)));
    end
    if (the1_Matrix1_less(i,j)<pi/2 && (phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))>pi/2)||((the1_Matrix1_less(i,j)<-pi/2 && (phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))>-pi/2))
        Ipeak_t2t3_Matrix1_less(i,j)=abs(Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)/Z1(i));

    else
        Ipeak_t2t3_Matrix1_less(i,j)=max(abs(Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*sin(the1_Matrix1_less(i,j))/Z1(i)),abs(Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*sin(phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))/Z1(i)));

    end
    Ipeak_Matrix1_less(i,j)=max(Ipeak_t0t2_Matrix1_less(i,j),Ipeak_t2t3_Matrix1_less(i,j));




if (abs(Vin_given(m(i,1))-Vin_Matrix1_less(i,j))<2)&&(s1==0)
    r0_Matrix_less(i,m(i,1))=x(1,i);
    phi0_Matrix_less(i,m(i,1))=x(2,i);
    the0_Matrix_less(i,m(i,1))=x(3,i);
    r1_Matrix_less(i,m(i,1))=x(4,i);
    phi1_Matrix_less(i,m(i,1))=x(5,i);
    the1_Matrix_less(i,m(i,1))=x(6,i);
    M_Matrix_less(i,m(i,1))=x(7,i); 
    Vin_Matrix_less(i,m(i,1))=Vo/n(i)/M_Matrix_less(i,m(i,1));
    fs_Matrix_less(i,m(i,1))=fs_Matrix1_less(i,j);

 %%%%%计算电流值%%%%%%%        ------包括输出电流平均值，谐振电流有效值以及谐振电流峰值

    ILs_rms_Matrix_less(i,m(i,1))=sqrt((r0_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*(2*phi0_Matrix1_less(i,j)+sin(2*the0_Matrix1_less(i,j))-sin(2*phi0_Matrix1_less(i,j)+2*the0_Matrix1_less(i,j)))/wr(i)/Ts(i,j)/Z0(i)/Z0(i)/2)+(r1_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*(2*phi1_Matrix1_less(i,j)+sin(2*the1_Matrix1_less(i,j))-sin(2*phi1_Matrix1_less(i,j)+2*the1_Matrix1_less(i,j)))/wn(i)/Ts(i,j)/Z1(i)/Z1(i)/2));
    if (the0_Matrix1_less(i,j)<pi/2 && (phi0_Matrix1_less(i,j)+the0_Matrix1_less(i,j))>pi/2)||((the0_Matrix1_less(i,j)<-pi/2 && (phi0_Matrix1_less(i,j)+the0_Matrix1_less(i,j))>-pi/2))
        Ipeak_t0t2_Matrix_less(i,m(i,1))=abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)/Z0(i));
    else
        Ipeak_t0t2_Matrix_less(i,m(i,1))=max(abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*sin(the0_Matrix1_less(i,j))/Z0(i)),abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*sin(phi0_Matrix1_less(i,j)+the0_Matrix1_less(i,j))/Z0(i)));
    end
    if (the1_Matrix1_less(i,j)<pi/2 && (phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))>pi/2)||((the1_Matrix1_less(i,j)<-pi/2 && (phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))>-pi/2))
        Ipeak_t2t3_Matrix_less(i,m(i,1))=abs(Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)/Z1(i));

    else
        Ipeak_t2t3_Matrix_less(i,m(i,1))=max(abs(Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*sin(the1_Matrix1_less(i,j))/Z1(i)),abs(Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*sin(phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))/Z1(i)));

    end
    Ipeak_Matrix_less(i,m(i,1))=max(Ipeak_t0t2_Matrix_less(i,m(i,1)),Ipeak_t2t3_Matrix_less(i,m(i,1)));
    
    m(i,1)=m(i,1)-1;
    if(m(i,1)==0)
        m(i,1)=1;
       s1=1;
    end
end

    %%%%%%   判断增益是否符合   %%%%%% 
    if M_Matrix1_less(i,j) > M_max(i)         %%% 欠谐振下的增益是大于1的
        Flag_M_less(i) = 1;
        aaa(i)=1;
        break
    end
    j=j+1;
end

%%%%%%      如果已经判定为不可行解，则直接跳过该循环，不再进行过谐振的迭代以及绘图
if Flag_error(i)==1 || Flag_M_less(i) == 0
    Flag_error(i)=1;
    continue
end


%%  过谐振--牛顿迭代法

%%%%%%  过谐振，不同频率范围表达式不一样    %%%%%%   

M(i)=1;
r0(i)=sqrt(Ir0N(i)^2+(Vcr0N(i)-(1-M(i)))^2);
the0(i)=atan(-Ir0N(i)/(Vcr0N(i)-(1-M(i))));
r1(i)=sqrt(Ir2N(i)^2+(Vcr2N(i)+1+M(i))^2);
the1(i)=pi+atan(-Ir2N(i)/(Vcr2N(i)+1+M(i)));
phi0(i)=t2(i)*wr0(i);
phi1(i)=tr1(i)*wr0(i);


%定义七个未知待求解参数r0,phi0,the0,r1,phi1,the1,M
x(1:7,i)=[r0(i),phi0(i),the0(i),r1(i),phi1(i),the1(i),M(i)]; % 7xPopsize

j=1;

for j_fs=(fr(i)):500:2.5*fr(i)
    if m(i,3)==11
        break%表示全工作在欠谐振区域内，不再进行过谐振的牛顿迭代求解
    end
    
    if Flag_error(i) == 1       %%%%%%      如果欠谐振出错，则直接不进行过谐振计算
        break
    end 

    %%%%%随着扫频改变周期%%%%%%%
    Ts(i,j)=1/j_fs;
    fs_Matrix1_more(i,j)=j_fs;
    iteration_times(i)=iteration_times(i)+1;
    
        for k=1:5000
        if k==5000
            Flag_newton(i)=1; %Flag_1=1代表收敛失败
            error_1(i) =100;
        break
        end

        b=F_more(x(:,i));%b输出是7个方程，方程理论上均=0

        norm_b=norm(b,inf);%矩阵的无穷范数是–各元素先取绝对值而后按行相加的最大值 ，本质上就是取最大的abs(A)，得到一个最大绝对值。
        if norm_b<eps1 %输出的7个方程中的最大值小于约定值则视为收敛完成跳出循环
            break;
        end
        A=Jac_more(x(:,i)); %A是导数
        dx=-pinv(A)*b; %见牛顿下山法公式(本文lamda没有下山逐渐减小，这样的优点是更快缺点是容易不收敛)
        x(:,i)=x(:,i)+lemda*dx; %新的迭代参数
        norm_dx=norm(dx,inf);%矩阵的无穷范数是–各元素先取绝对值而后按行相加的最大值
        if norm_dx<eps2
            break;
        end
        end

        %%%%%%  将迭代后的数据记录在矩阵当中      %%%%%%

        r0_Matrix1_more(i,j)=x(1,i);
        phi0_Matrix1_more(i,j)=x(2,i);
        the0_Matrix1_more(i,j)=x(3,i);
        r1_Matrix1_more(i,j)=x(4,i);
        phi1_Matrix1_more(i,j)=x(5,i);
        the1_Matrix1_more(i,j)=x(6,i);
        M_Matrix1_more(i,j)=x(7,i)    ;
        Vin_Matrix1_more(i,j)=Vo/n(i)/M_Matrix1_more(i,j);  %%%%%% 可能没必要
  
    
    %%%%%%  判断是否满足软开关条件     %%%%%%
            if abs(Vin_Matrix1_more(i,j)*r0_Matrix1_more(i,j)*sin(the0_Matrix1_more(i,j)+phi0_Matrix1_more(i,j))/Z0(i))<(2*Vin_Matrix1_more(i,j)*Coss/tdead)
%                 Flag_error(i)=1;
                Flag_ZVS(i)=1;
                error_1(i) =1;
                break
            end

        %%%%%%      判断是否进入NOP模态，即当前迭代结果是否可靠     %%%%%%
        
        Vlm_N(i,j)=Lm(i)/(Lm(i)+Lr(i))*(Vin_Matrix1_more(i,j)-Vin_Matrix1_more(i,j)*(-r0_Matrix1_more(i,j)*cos(the0_Matrix1_more(i,j))+1-M_Matrix1_more(i,j)));
        if Vlm_N <=Vo/n(i)
            Flag_NOP(i) =1;
            Flag_error(i)=1;
            error_1(i) =1;
            break
        end

        %%%%%计算电流值%%%%%%%        ------包括输出电流平均值，谐振Lr电流有效值以及谐振电流峰值
    Io_avg_Matrix1_more(i,j)=2*Vin_Matrix1_more(i,j)*(r0_Matrix1_more(i,j)*cos(the0_Matrix1_more(i,j))-r0_Matrix1_more(i,j)*cos(phi0_Matrix1_more(i,j)+the0_Matrix1_more(i,j))+r1_Matrix1_more(i,j)*cos(the1_Matrix1_more(i,j))-r1_Matrix1_more(i,j)*cos(phi1_Matrix1_more(i,j)+the1_Matrix1_more(i,j)))/n(i)/wr(i)/Ts(i,j)/Z0(i);
    ILs_rms_Matrix1_more(i,j)=sqrt((r0_Matrix1_more(i,j)*r0_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*(2*phi0_Matrix1_more(i,j)+sin(2*the0_Matrix1_more(i,j))-sin(2*phi0_Matrix1_more(i,j)+2*the0_Matrix1_more(i,j)))/wr(i)/Ts(i,j)/Z0(i)/Z0(i)/2)+(r1_Matrix1_more(i,j)*r1_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*(2*phi1_Matrix1_more(i,j)+sin(2*the1_Matrix1_more(i,j))-sin(2*phi1_Matrix1_more(i,j)+2*the1_Matrix1_more(i,j)))/wr(i)/Ts(i,j)/Z0(i)/Z0(i)/2));
    
    
    %%%%%%过谐振下的谐振电流最大值应该只出现在P模态的中间或者结尾处，因此判断相角是否过pi/2
    if (phi0_Matrix1_more(i,j)+the0_Matrix1_more(i,j))>pi/2 && (the0_Matrix1_more(i,j)<pi/2 )
        Ipeak_Matrix1_more(i,j)=abs(Vin_Matrix1_more(i,j)*r0_Matrix1_more(i,j)/Z0(i));
    else
        Ipeak_Matrix1_more(i,j)=abs(Vin_Matrix1_more(i,j)*r0_Matrix1_more(i,j)*sin(phi0_Matrix1_more(i,j)+the0_Matrix1_more(i,j))/Z0(i));
    end

if(abs(Vin_given(m(i,2)+1)-Vin_Matrix1_more(i,j))<2.5)&&(s2==0)
    r0_Matrix_more(i,m(i,2)+1)=x(1,i);
    phi0_Matrix_more(i,m(i,2)+1)=x(2,i);
    the0_Matrix_more(i,m(i,2)+1)=x(3,i);
    r1_Matrix_more(i,m(i,2)+1)=x(4,i);
    phi1_Matrix_more(i,m(i,2)+1)=x(5,i);
    the1_Matrix_more(i,m(i,2)+1)=x(6,i);
    M_Matrix_more(i,m(i,2)+1)=x(7,i); 
    Vin_Matrix_more(i,m(i,2)+1)=Vo/n(i)/M_Matrix_more(i,m(i,2)+1);
    fs_Matrix_more(i,m(i,2)+1)=fs_Matrix1_more(i,j);


%%%%%计算电流值%%%%%%%        ------包括输出电流平均值，谐振Lr电流有效值以及谐振电流峰值
    ILs_rms_Matrix_more(i,m(i,2)+1)=sqrt((r0_Matrix1_more(i,j)*r0_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*(2*phi0_Matrix1_more(i,j)+sin(2*the0_Matrix1_more(i,j))-sin(2*phi0_Matrix1_more(i,j)+2*the0_Matrix1_more(i,j)))/wr(i)/Ts(i,j)/Z0(i)/Z0(i)/2)+(r1_Matrix1_more(i,j)*r1_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*(2*phi1_Matrix1_more(i,j)+sin(2*the1_Matrix1_more(i,j))-sin(2*phi1_Matrix1_more(i,j)+2*the1_Matrix1_more(i,j)))/wr(i)/Ts(i,j)/Z0(i)/Z0(i)/2));
    
    
    %%%%%%过谐振下的谐振电流最大值应该只出现在P模态的中间或者结尾处，因此判断相角是否过pi/2
    if (phi0_Matrix1_more(i,j)+the0_Matrix1_more(i,j))>pi/2 && (the0_Matrix1_more(i,j)<pi/2 )
        Ipeak_Matrix_more(i,m(i,2)+1)=abs(Vin_Matrix1_more(i,j)*r0_Matrix1_more(i,j)/Z0(i));
    else
        Ipeak_Matrix_more(i,m(i,2)+1)=abs(Vin_Matrix1_more(i,j)*r0_Matrix1_more(i,j)*sin(phi0_Matrix1_more(i,j)+the0_Matrix1_more(i,j))/Z0(i));
    end

    m(i,2)=m(i,2)+1;
    if(m(i,2)==11)
          m(i,2)=10;
        s2=1;
    end
end

        %%%%%%   判断增益是否符合   %%%%%% 
        if M_Matrix1_more(i,j) < M_min(i)          %%% 过谐振下的增益是大于1的
        Flag_M_more(i) = 1;
        aaa(i)=1;
        break
        end
        j=j+1;
       
end

if Flag_M_more(i) == 0
    error_1(i) =1;
    Flag_error(i)=1;
end


%%%%%%      如果已经判定为不可行解，则直接跳过该循环，不再进行过谐振的迭代以及绘图
if Flag_error(i)==1 || Flag_M_more(i) == 0
    continue
end

 %%% 计算整个工作区域谐振电流有效值的最大值和最大峰值


%%%%%%  计算符合增益M条件且可行的频率范围   %%%%%%  ---默认从准谐振开始
if m(i,3)==11
    fs_min(i)=fs_Matrix_less(i,1);
    fs_max(i)=fs_Matrix_less(i,11);
    ILs_rms_max(i)=max(ILs_rms_Matrix1_less(i,:));
    Ipeak_max(i)=max(Ipeak_Matrix1_less(i,:));
    
else if m(i,3)==0
    fs_min(i)=fs_Matrix_more(i,1);
    fs_max(i)=fs_Matrix_more(i,11);
    ILs_rms_max(i)=max(ILs_rms_Matrix1_more(i,:));
    Ipeak_max(i)=max(Ipeak_Matrix1_more(i,:));
    else 
    fs_min(i)=fs_Matrix_less(i,1);
    fs_max(i)=fs_Matrix_more(i,11);
    ILs_rms_max(i)=max(max(ILs_rms_Matrix1_less(i,:)),max(ILs_rms_Matrix1_more(i,:)));
    Ipeak_max(i)=max(max(Ipeak_Matrix1_less(i,:)),max(Ipeak_Matrix1_more(i,:)));
    end
end



end


%% 欠谐振计算方程

function y=F_less(x) %%%%%通过七个代求解变量得到七个求解方程，且方程理论上都应为0%%%%%%%
% global Lr Cr Lm R Vin n Ts
% global Cr
y=zeros(size(x));
%%%%%待求解的7个变量%%%%%%%
r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);

%%%%%部分参数定义%%%%%%%
Ln(i)=Lm(i)/Lr(i);
Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i));wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
In(i)=Vin(i)/Z0(i);
IoN(i)=2*R*In(i)/(Ts(i,j)*wr0(i)*Vin(i)*n(i)*n(i));

%%%%%7个方程求解7个变量%%%%%%%
y(1)=r0*sin(phi0+the0)-r1/(sqrt(1+Ln(i)))*sin(the1);
y(2)=-r0*cos(phi0+the0)-M+r1*cos(the1);
y(3)=r1/(sqrt(1+Ln(i)))*sin(phi1+the1)+r0*sin(the0);
y(4)=-r1*cos(phi1+the1)-r0*cos(the0)+2-M;
y(5)=r0*sin(phi0+the0)-r0*sin(the0)-M/Ln(i)*phi0;
y(6)=M-IoN(i)*(r0*cos(the0)-r0*cos(phi0+the0)-r0*sin(the0)*phi0-M/(2*Ln(i))*phi0^2);
y(7)=phi0/wr0(i)+phi1/wr1(i)-Ts(i,j)/2;

end

%%  欠谐振方程导数

function Mx=Jac_less(x)

%%%%%待求解的7个变量%%%%%%%
r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);

%%%%%部分参数定义%%%%%%%

Ln(i)=Lm(i)/Lr(i);
Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i));wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
In(i)=Vin(i)/Z0(i);
IoN(i)=2*R*In(i)/(Ts(i,j)*wr0(i)*Vin(i)*n(i)*n(i));

%%%%%对方程针对每个待求参数求偏导%%%%%%%
Mx=[sin(phi0+the0) r0*cos(phi0+the0) r0*cos(phi0+the0) -1/sqrt(1+Ln(i))*sin(the1) 0 r1/sqrt(1+Ln(i))*cos(the1) 0;...
    -cos(phi0+the0) r0*sin(phi0+the0) r0*sin(phi0+the0) cos(the1) 0 -r1*sin(the1) -1;...
    sin(the0) 0 r0*cos(the0) 1/sqrt(1+Ln(i))*sin(phi1+the1) r1/sqrt(1+Ln(i))*cos(phi1+the1) r1/sqrt(1+Ln(i))*cos(phi1+the1) 0;...
    -cos(the0) 0 r0*sin(the0) -cos(phi1+the1) r1*sin(phi1+the1) r1*sin(phi1+the1) -1;...
    sin(phi0+the0)-sin(the0) r0*cos(phi0+the0)-M/Ln(i) r0*cos(phi0+the0)-r0*cos(the0) 0 0 0 -phi0/Ln(i);...
    -IoN(i)*(cos(the0)-cos(phi0+the0)-sin(the0)*phi0) -IoN(i)*(r0*sin(phi0+the0)-r0*sin(the0)-M/Ln(i)*phi0) -IoN(i)*(-r0*sin(the0)+r0*sin(phi0+the0)-r0*cos(the0)*phi0) 0 0 0 1-IoN(i)*(-1/(2*Ln(i))*phi0^2);...
    0 1/wr0(i) 0 0 1/wr1(i) 0 0]; 


end

%% 过谐振计算方程

function y=F_more(x) %%%%%通过七个代求解变量得到七个求解方程，且方程理论上都应为0%%%%%%%
% global Lr Cr Lm R Vin n Ts i
% global Cr
y=zeros(size(x));
%%%%%待求解的7个变量%%%%%%%
r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);

%%%%%部分参数定义%%%%%%%
% Lr=22e-6;Cr=114e-9;Lm=80e-6;R=40;Vin=60;n=1.04; Ts=13.3333e-6;
Ln(i)=Lm(i)/Lr(i);
Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i));wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
In(i)=Vin(i)/Z0(i);
IoN(i)=2*R*In(i)/(Ts(i,j)*wr0(i)*Vin(i)*n(i)*n(i));

%%%%%7个方程求解7个变量%%%%%%%
y(1)=r0*sin(the0)+M*wr0(i)*Ts(i,j)/(4*Ln(i));
y(2)=r0*sin(phi0+the0)-r1*sin(the1);
y(3)=-r0*cos(phi0+the0)+r1*cos(the1)+2;
y(4)=r1*sin(phi1+the1)+r0*sin(the0);
y(5)=-r1*cos(phi1+the1)-r0*cos(the0)-2*M;
y(6)=M-IoN(i)*(r0*cos(the0)-r0*cos(phi0+the0)+r1*cos(the1)-r1*cos(phi1+the1));
y(7)=phi0/wr0(i)+phi1/wr0(i)-Ts(i,j)/2;


end


%% 过谐振的方程导数
function Mx=Jac_more(x)

% global Lr Cr Lm R Vin n Ts i
% global Cr
%%%%%待求解的7个变量%%%%%%%
r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);


%%%%%部分参数定义%%%%%%%
% Lr=22e-6;Cr=114e-9;Lm=80e-6;R=40;Vin=60;n=1.04; Ts=13.3333e-6;
Ln(i)=Lm(i)/Lr(i);
Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i));wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
In(i)=Vin(i)/Z0(i);
IoN(i)=2*R*In(i)/(Ts(i,j)*wr0(i)*Vin(i)*n(i)*n(i));

%%%%%对方程针对每个待求参数求偏导%%%%%%%



Mx=[sin(the0) 0 r0*cos(the0) 0 0 0 Ts(i,j)*wr0(i)/(4*Ln(i));...
    sin(phi0+the0) r0*cos(phi0+the0) r0*cos(phi0+the0) -sin(the1) 0 -r1*cos(the1) 0;...
    -cos(phi0+the0) r0*sin(phi0+the0) r0*sin(phi0+the0) cos(the1) 0 -r1*sin(the1) 0;...
    sin(the0) 0 r0*cos(the0) sin(phi1+the1) r1*cos(phi1+the1) r1*cos(phi1+the1) 0;...
    -cos(the0) 0 r0*sin(the0) -cos(phi1+the1) r1*sin(phi1+the1) r1*sin(phi1+the1) -2;...
    -IoN(i)*(cos(the0)-cos(phi0+the0)) -IoN(i)*r0*sin(phi0+the0) -IoN(i)*(-r0*sin(the0)+r0*sin(phi0+the0)) -IoN(i)*(cos(the1)-cos(phi1+the1)) -IoN(i)*r1*sin(phi1+the1) -IoN(i)*(-r1*sin(the0)+r1*sin(phi1+the1)) 1;...
    0 1/wr0(i) 0 0 1/wr0(i) 0 0]; 

end



1
end