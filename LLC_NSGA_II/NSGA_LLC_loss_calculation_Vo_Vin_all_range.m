function[P_given Pswitch_off_avg,P_conduction_Mosfet_avg,P_conduction_Diode_avg,P_Tr_Cu_avg,P_Lr_Cu_avg,P_Tr_Fe_avg,P_Lr_Fe_avg, P_recovery_Diode_Matrix_avg,Ploss_Cell ,Pswitch_off_Cell,P_conduction_Mosfet_Cell,P_conduction_Diode_Cell,P_Tr_Cu_Cell,P_Lr_Cu_Cell,Bm_T_t0t2_Cell,P_Tr_Fe_t0t2_Cell,P_Tr_Fe_t2t3_Cell,P_Lr_Fe_t0t2_Cell,P_Lr_Fe_t2t3_Cell,ILs_rms_Cell,Im1_Cell,fs_Cell,Vin_Cell,Vin_Matrix1_Cell ] = NSGA_LLC_loss_calculation_Vo_Vin_all_range(n,Lr,fr,Lm,R,Rdc1,Rdc2,Rdc3,N1,N2,N3,Ae_T,Ae_L,Ve_T,Ve_L,Flag_error)


global Vo Vin_max Vin_min
global Vin_given fs_min fs_max Cr Ln Z0 Z1 wr wn ILs_rms_max Ipeak_max ...
ILs_rms_Matrix_less Ipeak_t0t2_Matrix_less Ipeak_t2t3_Matrix_less Ipeak_Matrix_less fs_Matrix_less  ...
r0_Matrix_less phi0_Matrix_less the0_Matrix_less  ...
r1_Matrix_less phi1_Matrix_less the1_Matrix_less M_Matrix_less Vin_Matrix_less  ...
ILs_rms_Matrix_more Ipeak_Matrix_more fs_Matrix_more  r0_Matrix_more  phi0_Matrix_more the0_Matrix_more  ...
r1_Matrix_more phi1_Matrix_more the1_Matrix_more M_Matrix_more Vin_Matrix_more  

global M_max  M_min 



global pop_size tf Ron alpha beta kcof Vo Vin_max Izvs  Qrr

%%%%%%%%%%%%%参数初始化%%%%%%%%%%%%%%


%%%%%%%%%%%%%   针对欠谐振的参数变量的定义及初始化   %%%%%%%%%%%

Im1_Matrix_less=zeros(pop_size,11);
Pswitch_off_Matrix_less=zeros(pop_size,11);
P_conduction_Mosfet_Matrix_less=zeros(pop_size,11);
P_conduction_Diode_Matrix_less=zeros(pop_size,11);
P_Tr_Cu_Matrix_less=zeros(pop_size,11);
P_Lr_Cu_Matrix_less=zeros(pop_size,11);
Bm_T_t0t2_less=zeros(pop_size,11);
Bm_T_t2t3_less=zeros(pop_size,11);
P_Tr_Fe_t0t2_Matrix_less=zeros(pop_size,11);
P_Tr_Fe_t2t3_Matrix_less=zeros(pop_size,11);
Bm_L_t0t2_less=zeros(pop_size,11);
P_Lr_Fe_t0t2_Matrix_less=zeros(pop_size,11);
P_Lr_Fe_t2t3_Matrix_less=zeros(pop_size,11);
Bm_L_t2t3_less=zeros(pop_size,11);



%%%%%%%%%%%%%   针对过谐振的参数变量的定义及初始化   %%%%%%%%%%%

Vf=0.9;


Im1_Matrix_more=zeros(pop_size,11);
Pswitch_off_Matrix_more=zeros(pop_size,11);
P_conduction_Mosfet_Matrix_more=zeros(pop_size,11);
P_conduction_Diode_Matrix_more=zeros(pop_size,11);
P_Tr_Cu_Matrix_more=zeros(pop_size,11);
P_Lr_Cu_Matrix_more=zeros(pop_size,11);

Bm_T_more=zeros(pop_size,11);
P_Tr_Fe_Matrix_more=zeros(pop_size,11);
Bm_L_t0t2_more=zeros(pop_size,11);
P_Lr_Fe_t0t2_Matrix_more=zeros(pop_size,11);
P_Lr_Fe_t2t3_Matrix_more=zeros(pop_size,11);
Bm_L_t2t3_more=zeros(pop_size,11);
P_recovery_Diode_Matrix_more=zeros(pop_size,11);


%%%%%%%         损耗平均值定义         %%%%%%%%%



Pswitch_off_avg=zeros(pop_size,1);
P_conduction_Mosfet_avg=zeros(pop_size,1);
P_conduction_Diode_avg=zeros(pop_size,1);
P_recovery_Diode_Matrix_avg=zeros(pop_size,1);
P_Tr_Cu_avg=zeros(pop_size,1);
P_Lr_Cu_avg=zeros(pop_size,1);
P_Tr_Fe_t0t2_avg=zeros(pop_size,1);
P_Tr_Fe_t2t3_avg=zeros(pop_size,1);
P_Tr_Fe_avg=zeros(pop_size,1);
P_Lr_Fe_t0t2_avg=zeros(pop_size,1);
P_Lr_Fe_t2t3_avg=zeros(pop_size,1);
P_Lr_Fe_avg=zeros(pop_size,1);



Ploss_Cell=cell(1,pop_size);
Pswitch_off_Cell=cell(1,pop_size);
P_conduction_Mosfet_Cell=cell(1,pop_size);
P_conduction_Diode_Cell=cell(1,pop_size);
Bm_T_t0t2_Cell=cell(1,pop_size);
P_Tr_Cu_Cell=cell(1,pop_size);
P_Lr_Cu_Cell=cell(1,pop_size);
P_Tr_Fe_t0t2_Cell=cell(1,pop_size);
P_Tr_Fe_t2t3_Cell=cell(1,pop_size);
P_Lr_Fe_t0t2_Cell=cell(1,pop_size);
P_Lr_Fe_t2t3_Cell=cell(1,pop_size);
ILs_rms_Cell=cell(1,pop_size);
Im1_Cell=cell(1,pop_size);
fs_Cell=cell(1,pop_size);
Vin_Cell=cell(1,pop_size);
Vin_Matrix1_Cell=cell(1,pop_size);


%%%%%%% 定义特定输入电压下的效率
P_given=zeros(pop_size,11);

factor=2;

    %%%%%%   损耗迭代计算 %%%%%%  
    for i=1:pop_size
        if n(i)<=Vo/Vin_max
            m(i)=11;
        else if n(i)>=Vo/Vin_min
                m(i)=0;
            
        else                 
        for j=1:1:10
            if ((n(i)*Vin_given(j))<=Vo)&&((n(i)*Vin_given(j+1))>Vo)
                m(i)=j;
                break
            end
        end
        end      
        end
    end


for i=1:pop_size
    

    if Flag_error(i)==0

            
            %%%%%%  欠谐振迭代计算 %%%%%%%
            if m(i)~=0
                
           


            %%%%%%%%%%%%%损耗计算部分%%%%%%%%%%%%%%
            for j=1:m(i)  
                %%%%%%%%%%%%%%开关损耗%%%%%%%%%%%%%%       

                Im1_Matrix_less(i,j)=Vin_Matrix_less(i,j)*r0_Matrix_less(i,j)*sin(the0_Matrix_less(i,j))/Z0(i);
               
                Pswitch_off_Matrix_less(i,j)=4*Vin_Matrix_less(i,j)*abs(Im1_Matrix_less(i,j))*tf*fs_Matrix_less(i,j)/2; %原边开关管关断损耗 
                
                %%%%%%%%%%%%%%导通损耗%%%%%%%%%%%%%%

                %——MOSFET导通损耗——%

                P_conduction_Mosfet_Matrix_less(i,j)=2*Ron*ILs_rms_Matrix_less(i,j)*ILs_rms_Matrix_less(i,j);%MOSFET导通损耗 一个周期内同时存在两个导通
                %——二极管导通损耗——%
                            
                % Conduction = @(phi) 0.7 *power(Vin_Matrix_less(i,j)*(r0_Matrix_less(i,j)*sin(phi+the0_Matrix_less(i,j))-r0_Matrix_less(i,j)*sin(the0_Matrix_less(i,j))-M_Matrix_less(i,j)*Lr(i)*phi/Lm(i))/Z0(i)/n(i),0.1861444575812+1);
                % P_conduction_Diode_Matrix_less(i,j)=1/phi0_Matrix_less(i,j)*2*integral(Conduction,0,phi0_Matrix_less(i,j));  %二极管导通损耗 一个周期内同时存在两个导通 
                P_conduction_Diode_Matrix_less(i,j) = Vf * 5 * 2;

                %%%%%%%%%%%%%% 计算交流等效电阻
                Rac1(i,j)=Rdc1(i)*(1.5+(fs_Matrix_less(i,j)-20e3)/(350e3-20e3)*2.5);
                Rac2(i,j)=Rdc2(i)*(1.5+(fs_Matrix_less(i,j)-20e3)/(350e3-20e3)*2.5);
                Rac3(i,j)=Rdc3(i)*(1.5+(fs_Matrix_less(i,j)-20e3)/(350e3-20e3)*2.5);
%                 Rac1(i,j)=Rdc1(i);
%                 Rac2(i,j)=Rdc2(i);
%                 Rac3(i,j)=Rdc3(i);

                %%%%%%%%%%%%%%  铜损  %%%%%%%%%%%%%%
                P_Tr_Cu_Matrix_less(i,j)=ILs_rms_Matrix_less(i,j)*ILs_rms_Matrix_less(i,j)*Rac1(i,j)...
                    +ILs_rms_Matrix_less(i,j)*ILs_rms_Matrix_less(i,j)*Rac2(i,j)/n(i)/n(i);    %变压器铜损
                P_Lr_Cu_Matrix_less(i,j)=ILs_rms_Matrix_less(i,j)*ILs_rms_Matrix_less(i,j)*Rac3(i,j);                         %谐振电感铜损
               
                %%%%%%%%%%%%%%  磁损  %%%%%%%%%%%%%%

                  %%t0t2
                ir_t0_less(i,j)=Vin_Matrix_less(i,j)*r0_Matrix_less(i,j)*sin(the0_Matrix_less(i,j))/Z0(i);
                ir_t2_less(i,j)=Vin_Matrix_less(i,j)*r0_Matrix_less(i,j)*sin(the0_Matrix_less(i,j)+phi0_Matrix_less(i,j))/Z0(i);
                deta_B_T_t0t2_less(i,j)=Lm(i)*(ir_t2_less(i,j)-ir_t0_less(i,j))/N1(i)/Ae_T(i);
                
                dB_T_t0t2_less(i,j)= Vo/N2(i)/Ae_T(i);% P模态，电压励磁，磁密的变化率一定，其计算公式可以等效为变压器二次侧计算
            
                P_Tr_Fe_t0t2_Matrix_less(i,j)=(2*fs_Matrix_less(i,j)*phi0_Matrix_less(i,j)/wr(i))*kcof*Ve_T(i)... 
                    *power(dB_T_t0t2_less(i,j),alpha)*power( deta_B_T_t0t2_less(i,j),(beta-alpha));%系数2表示开关周期当中有两个P模态

                %%t2t3
   
                ir_t3_less(i,j)=Vin_Matrix_less(i,j)*r1_Matrix_less(i,j)*sin(the1_Matrix_less(i,j)+phi1_Matrix_less(i,j))/Z1(i);
            
                deta_i_t2t3_less(i,j)=max(abs(Ipeak_t2t3_Matrix_less(i,j)-ir_t2_less(i,j)),abs(Ipeak_t2t3_Matrix_less(i,j)-ir_t3_less(i,j)));
                deta_B_T_t2t3_less(i,j)=Lm(i)*(deta_i_t2t3_less(i,j))/N1(i)/Ae_T(i);
    
                fun1_less = @(phi) power(abs(Lm(i)*Vin_Matrix_less(i,j)*r1_Matrix_less(i,j)*wn(i)*cos(phi+the1_Matrix_less(i,j))/Z1(i)/N1(i)/Ae_T(i)),alpha);
                
                P_Tr_Fe_t2t3_Matrix_less(i,j)=2*fs_Matrix_less(i,j)...
                *power(deta_B_T_t2t3_less(i,j),beta-alpha)*integral(fun1_less,0,phi1_Matrix_less(i,j))*kcof*Ve_T(i)/wn(i);

                %——PO模态谐振电感磁损——%
                deta_B_L_t0t2_less(i,j)=Lr(i)*(Ipeak_t0t2_Matrix_less(i,j)-ir_t0_less(i,j))/N3(i)/Ae_L(i);
                deta_B_L_t2t3_less(i,j)=Lr(i)*(deta_i_t2t3_less(i,j))/N3(i)/Ae_L(i);;


                fun2_less = @(phi) power(abs(Lr(i)*Vin_Matrix_less(i,j)*r0_Matrix_less(i,j)*wr(i) ...
                    *cos(phi+the0_Matrix_less(i,j))/N3(i)/Ae_L(i)/Z0(i)),alpha);

                P_Lr_Fe_t0t2_Matrix_less(i,j)=2*fs_Matrix_less(i,j)*power(deta_B_L_t0t2_less(i,j),beta-alpha) ...
                    *integral(fun2_less,0,phi0_Matrix_less(i,j))*kcof*Ve_L(i)/wr(i);

                fun3_less = @(phi) power(abs(Lr(i)*Vin_Matrix_less(i,j)*r1_Matrix_less(i,j)*wn(i) ...
                    *cos(phi+the1_Matrix_less(i,j))/N3(i)/Ae_L(i)/Z1(i)),alpha);

                P_Lr_Fe_t2t3_Matrix_less(i,j)=2*fs_Matrix_less(i,j)*power( deta_B_L_t2t3_less(i,j),beta-alpha) ...
                    *integral(fun3_less,0,phi1_Matrix_less(i,j))*kcof*Ve_L(i)/wn(i);

                
                P_given(i,j)=Pswitch_off_Matrix_less(i,j)+P_conduction_Mosfet_Matrix_less(i,j)+P_conduction_Diode_Matrix_less(i,j)+P_Tr_Cu_Matrix_less(i,j)+P_Lr_Cu_Matrix_less(i,j)...
                    +P_Tr_Fe_t0t2_Matrix_less(i,j)+P_Tr_Fe_t2t3_Matrix_less(i,j)+P_Lr_Fe_t0t2_Matrix_less(i,j)+P_Lr_Fe_t2t3_Matrix_less(i,j);

            end
            end

            %%       %%%%%%      过谐振损耗计算         %%%%%%

                %%%%%%%%%%%%%参数初始化%%%%%%%%%%%%%%
        
            if m(i)~= 11
            
               %%%%%%%%%%%%%损耗计算部分%%%%%%%%%%%%%%

            for j=(m(i)+1):11  
                %%%%%%%%%%%%%%开关损耗%%%%%%%%%%%%%%       

                Im1_Matrix_more(i,j)=Vin_Matrix_more(i,j)*r0_Matrix_more(i,j)*sin(the0_Matrix_more(i,j)+phi0_Matrix_more(i,j))/Z0(i);
               
                Pswitch_off_Matrix_more(i,j)=4*Vin_Matrix_more(i,j)*abs(Im1_Matrix_more(i,j))*tf*fs_Matrix_more(i,j)/2; %原边开关管关断损耗 
                

                %%%%%%%%%%%%%%导通损耗%%%%%%%%%%%%%%

                %——MOSFET导通损耗——%

                P_conduction_Mosfet_Matrix_more(i,j)=2*Ron*ILs_rms_Matrix_more(i,j)*ILs_rms_Matrix_more(i,j);%MOSFET导通损耗 一个周期内同时存在两个导通
               
                %——二极管导通损耗——%
                            
                %%%%%%  下面的计算是电流的0.2233次幂乘以系数0.56，然后积分

                %%%对fai求导？不等于对时间求导，相差角速度系数，，，2pi
                
                % Conduction_more_t0t2 = @(phi) 0.7*power(Vin_Matrix_more(i,j)*(r0_Matrix_more(i,j)*sin(phi+the0_Matrix_more(i,j))-r0_Matrix_more(i,j)*sin(the0_Matrix_more(i,j))-M_Matrix_more(i,j)*Lr(i)*phi/Lm(i))/Z0(i)/n(i),0.1861444575812+1);
                % Conduction_more_t2t3 = @(phi) 0.7*power(Vin_Matrix_more(i,j)*(r1_Matrix_more(i,j)*sin(phi+the1_Matrix_more(i,j))-r0_Matrix_more(i,j)*sin(the0_Matrix_more(i,j))-M_Matrix_more(i,j)*Lr(i)*(phi+phi0_Matrix_more(i,j))/Lm(i))/Z0(i)/n(i),0.1861444575812+1);
                % P_conduction_Diode_Matrix_more(i,j)=1/phi0_Matrix_more(i,j)*(2*integral(Conduction_more_t0t2,0,phi0_Matrix_more(i,j)))+1/phi1_Matrix_more(i,j)*(2*integral(Conduction_more_t2t3,0,phi1_Matrix_more(i,j)));  %二极管导通损耗 一个周期内同时存在两个导通 
                
                P_conduction_Diode_Matrix_more(i,j) = Vf * 5 *2;
               
                %%%%%%%%%%%%%%二极管反向恢复损耗%%%%%%%%%%%%%%
                
                P_recovery_Diode_Matrix_more(i,j)=4*Vo*Qrr*fs_Matrix_more(i,j);         %%%%%%%%%   一个周期内有四个二极管，采用全桥整流

                %%%%%%%%%%%%%% 计算交流等效电阻
                Rac1(i,j)=Rdc1(i)*(1.5+(fs_Matrix_more(i,j)-20e3)/(350e3-20e3)*2.5);
                Rac2(i,j)=Rdc2(i)*(1.5+(fs_Matrix_more(i,j)-20e3)/(350e3-20e3)*2.5);
                Rac3(i,j)=Rdc3(i)*(1.5+(fs_Matrix_more(i,j)-20e3)/(350e3-20e3)*2.5);

%                     Rac1(i,j)=Rdc1(i);
%                     Rac2(i,j)=Rdc2(i);
%                     Rac3(i,j)=Rdc3(i);


                %%%%%%%%%%%%%%  铜损  %%%%%%%%%%%%%%
                P_Tr_Cu_Matrix_more(i,j)=ILs_rms_Matrix_more(i,j)*ILs_rms_Matrix_more(i,j)*Rac1(i,j)+ILs_rms_Matrix_more(i,j)*ILs_rms_Matrix_more(i,j)*Rac2(i,j)/n(i)/n(i);    %变压器铜损
                
                P_Lr_Cu_Matrix_more(i,j)=ILs_rms_Matrix_more(i,j)*ILs_rms_Matrix_more(i,j)*Rac3(i,j);                         %谐振电感铜损
               
                %%%%%%%%%%%%%%  磁损  %%%%%%%%%%%%%%

                
                 %——NP变压器磁损——%
                %%%%过谐振NP模态下变压器都是电压励磁，可以不用区分N模态还是P模态
 
                deta_B_T_more(i,j)= Vo/N2(i)/2/fs_Matrix_more(i,j)/Ae_T(i);       
                dB_T_more(i,j)= Vo/N2(i)/Ae_T(i);% P/N模态，电压励磁，磁密的变化率一定，其计算公式可以等效为变压器二次侧计算
                P_Tr_Fe_Matrix_more(i,j)=(2*fs_Matrix_more(i,j)*pi/wr(i))*((power(dB_T_more(i,j),alpha)*power(deta_B_T_more(i,j),beta-alpha))*kcof*Ve_T(i));             

                %——NP谐振电感磁损——%
                ir_t0_more(i,j)=Vin_Matrix_more(i,j)*r0_Matrix_more(i,j)*sin(the0_Matrix_more(i,j))/Z0(i);
                ir_t2_more(i,j)=Vin_Matrix_more(i,j)*r0_Matrix_more(i,j)*sin(phi0_Matrix_more(i,j)+the0_Matrix_more(i,j))/Z0(i);
                ir_t3_more(i,j)=Vin_Matrix_more(i,j)*r1_Matrix_more(i,j)*sin(phi1_Matrix_more(i,j)+the1_Matrix_more(i,j))/Z0(i);

                deta_B_L_t0t2_more(i,j)=Lr(i)*( Ipeak_Matrix_more(i,j)-ir_t0_more(i,j))/N3(i)/Ae_L(i);
                deta_B_L_t2t3_more(i,j)=Lr(i)*(ir_t2_more(i,j)-ir_t3_more(i,j))/N3(i)/Ae_L(i);

                fun2_more = @(phi) power(abs(Lr(i)*Vin_Matrix_more(i,j)*r0_Matrix_more(i,j)*wr(i) ...
                    *cos(phi+the0_Matrix_more(i,j))/N3(i)/Ae_L(i)/Z0(i)),alpha);
                P_Lr_Fe_t0t2_Matrix_more(i,j)=2*fs_Matrix_more(i,j)*power( deta_B_L_t0t2_more(i,j),beta-alpha) ...
                    *integral(fun2_more,0,phi0_Matrix_more(i,j))*kcof*Ve_L(i)/wr(i);

                fun3_more = @(phi) power(abs(Lr(i)*Vin_Matrix_more(i,j)*r1_Matrix_more(i,j)*wr(i) ...
                    *cos(phi+the1_Matrix_more(i,j))/N3(i)/Ae_L(i)/Z0(i)),alpha);
                P_Lr_Fe_t2t3_Matrix_more(i,j)=2*fs_Matrix_more(i,j)*power(deta_B_L_t2t3_more(i,j),beta-alpha) ...
                    *integral(fun3_more,0,phi1_Matrix_more(i,j))*kcof*Ve_L(i)/wr(i);

                P_given(i,j)=Pswitch_off_Matrix_more(i,j)+P_conduction_Mosfet_Matrix_more(i,j)+ P_conduction_Diode_Matrix_more(i,j)+P_recovery_Diode_Matrix_more(i,j)+P_Tr_Cu_Matrix_more(i,j)+P_Lr_Cu_Matrix_more(i,j)+ P_Tr_Fe_Matrix_more(i,j)+P_Lr_Fe_t0t2_Matrix_more(i,j)+ P_Lr_Fe_t2t3_Matrix_more(i,j);  
            end
            else
                m(i)=11;
            end


 
%             Ploss_Cell(i)={[Pswitch_off_Matrix_less(i,1:m(i) ) Pswitch_off_Matrix_more(i,(m(i)+1):11)]+[P_conduction_Mosfet_Matrix_less(i,1:m(i)) P_conduction_Mosfet_Matrix_more(i,(m(i)+1):11)]+ [P_conduction_Diode_Matrix_less(i,1:m(i)) P_conduction_Diode_Matrix_more(i,(m(i)+1):11)]+[P_Tr_Cu_Matrix_less(i,1:m(i)) P_Tr_Cu_Matrix_more(i,(m(i)+1):11)]+[P_Lr_Cu_Matrix_less(i,1:m(i)) P_Lr_Cu_Matrix_more(i,(m(i)+1):11)]+[P_Tr_Fe_t0t2_Matrix_less(i,1:m(i))+P_Tr_Fe_t2t3_Matrix_less(i,m(i)) P_Tr_Fe_Matrix_more(i,(m(i)+1):11)]...
%                 +[P_Lr_Fe_t0t2_Matrix_less(i,1:m(i)) P_Lr_Fe_t0t2_Matrix_more(i,(m(i)+1):11)]+[P_Lr_Fe_t2t3_Matrix_less(i,1:m(i)) P_Lr_Fe_t2t3_Matrix_more(i,(m(i)+1):11)]};
            Pswitch_off_Cell(i)={[Pswitch_off_Matrix_less(i,1:m(i)) Pswitch_off_Matrix_more(i,(m(i)+1):11)]};
            P_conduction_Mosfet_Cell(i)={[P_conduction_Mosfet_Matrix_less(i,1:m(i)) P_conduction_Mosfet_Matrix_more(i,(m(i)+1):11)]};
            P_conduction_Diode_Cell(i)={[P_conduction_Diode_Matrix_less(i,1:m(i)) P_conduction_Diode_Matrix_more(i,(m(i)+1):11)]};
            P_Tr_Cu_Cell(i)={[P_Tr_Cu_Matrix_less(i,1:m(i)) P_Tr_Cu_Matrix_more(i,(m(i)+1):11)]};  
            P_Lr_Cu_Cell(i)={[P_Lr_Cu_Matrix_less(i,1:m(i)) P_Lr_Cu_Matrix_more(i,(m(i)+1):11)]};
            Bm_T_t0t2_Cell(i)={Bm_T_t0t2_less(i,1:m(i))};
            P_Tr_Fe_t0t2_Cell(i)={P_Tr_Fe_t0t2_Matrix_less(i,1:m(i)) };
            P_Tr_Fe_t2t3_Cell(i)={P_Tr_Fe_t2t3_Matrix_less(i,1:m(i))};
            P_Lr_Fe_t0t2_Cell(i)={[P_Lr_Fe_t0t2_Matrix_less(i,1:m(i)) P_Lr_Fe_t0t2_Matrix_more(i,(m(i)+1):11)]};
            P_Lr_Fe_t2t3_Cell(i)={[P_Lr_Fe_t2t3_Matrix_less(i,1:m(i)) P_Lr_Fe_t2t3_Matrix_more(i,(m(i)+1):11)]};

            ILs_rms_Cell(i)={[ILs_rms_Matrix_less(i,1:m(i)) ILs_rms_Matrix_more(i,(m(i)+1):11)]}; 

            Im1_Cell(i)={[Im1_Matrix_less(i,1:m(i)) Im1_Matrix_more(i,(m(i)+1):11)]}; 
            fs_Cell(i)={[fs_Matrix_less(i,1:m(i)) fs_Matrix_more(i,(m(i)+1):11)]}; 
            Vin_Cell(i)={[Vin_Matrix_less(i,1:m(i)) Vin_Matrix_more(i,(m(i)+1):11)] }; 
            Vin_Matrix1_Cell(i)={[Vin_Matrix_less(i,1:m(i)) Vin_Matrix_more(i,(m(i)+1):11)]}; 


            %%%%%%%%%%%%%-----平均损耗计算部分-----%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%开关损耗%%%%%%%%%%%%%%
            Pswitch_off_avg(i)=(sum(Pswitch_off_Matrix_less(i,1:m(i)))+sum(Pswitch_off_Matrix_more(i,(m(i)+1):11)))/11;  %原边开关管关断损耗 
            %%%%%%%%%%%%%%导通损耗%%%%%%%%%%%%%%
            P_conduction_Mosfet_avg(i)=(sum(P_conduction_Mosfet_Matrix_less(i,1:m(i)))+sum(P_conduction_Mosfet_Matrix_more(i,(m(i)+1):11)))/11 ;%MOSFET导通损耗 一个周期内同时存在两个导通
            P_conduction_Diode_avg(i)=(sum(P_conduction_Diode_Matrix_less(i,1:m(i)))+sum(P_conduction_Diode_Matrix_more(i,(m(i)+1):11)))/11;
            
            %%%%%%%%%%%%%%二极管反向恢复损耗%%%%%%%%%%%%%%

            P_recovery_Diode_Matrix_avg(i)=sum(P_recovery_Diode_Matrix_more(i,(m(i)+1):11))/11;

            %%%%%%%%%%%%%%  铜损  %%%%%%%%%%%%%%
            P_Tr_Cu_avg(i)=(sum(P_Tr_Cu_Matrix_less(i,1:m(i)))+sum(P_Tr_Cu_Matrix_more(i,(m(i)+1):11)))/11;  %变压器铜损
            P_Lr_Cu_avg(i)=(sum(P_Lr_Cu_Matrix_less(i,1:m(i)))+sum(P_Lr_Cu_Matrix_more(i,(m(i)+1):11)))/11;  %谐振电感铜损
            
            %%%%%%%%%%%%%%  磁损  %%%%%%%%%%%%%%
            
            %——变压器磁损——%
            
            P_Tr_Fe_avg(i)=(sum(P_Tr_Fe_t0t2_Matrix_less(i,1:m(i)))+sum(P_Tr_Fe_t2t3_Matrix_less(i,1:m(i)))+sum(P_Tr_Fe_Matrix_more(i,(m(i)+1):11)))/11;
      
            %——谐振电感磁损——%
           P_Lr_Fe_t0t2_avg(i)=(sum(P_Lr_Fe_t0t2_Matrix_less(i,1:m(i)))+sum(P_Lr_Fe_t0t2_Matrix_more(i,(m(i)+1):11)))/11;            
           P_Lr_Fe_t2t3_avg(i)=(sum(P_Lr_Fe_t2t3_Matrix_less(i,1:m(i)))+sum(P_Lr_Fe_t2t3_Matrix_more(i,(m(i)+1):11)))/11;
                  
            P_Lr_Fe_avg(i)=P_Lr_Fe_t0t2_avg(i)+P_Lr_Fe_t2t3_avg(i);
    end
end

%  Pswitch_off_avg
%  P_conduction_Mosfet_avg
%   P_conduction_Diode_avg
%   P_Tr_Cu_avg
%   P_Lr_Cu_avg
%   P_Tr_Fe_avg
%   P_Lr_Fe_avg

%%%%%%%%%%%%%-----转换参数存储为cell-----%%%%%%%%%%%%%%
Ploss_Cell=Ploss_Cell';
Pswitch_off_Cell=Pswitch_off_Cell';
P_conduction_Mosfet_Cell=P_conduction_Mosfet_Cell';
P_conduction_Diode_Cell=P_conduction_Diode_Cell';
P_Tr_Cu_Cell=P_Tr_Cu_Cell';
P_Lr_Cu_Cell=P_Lr_Cu_Cell';
Bm_T_t0t2_Cell=Bm_T_t0t2_Cell';
P_Tr_Fe_t0t2_Cell=P_Tr_Fe_t0t2_Cell';
P_Tr_Fe_t2t3_Cell=P_Tr_Fe_t2t3_Cell';
P_Lr_Fe_t0t2_Cell=P_Lr_Fe_t0t2_Cell';
P_Lr_Fe_t2t3_Cell=P_Lr_Fe_t2t3_Cell';
ILs_rms_Cell=ILs_rms_Cell';
Im1_Cell=Im1_Cell';
fs_Cell=fs_Cell';
Vin_Cell=Vin_Cell';
Vin_Matrix1_Cell=Vin_Matrix1_Cell';
1

% %仅用于对最后Pareto前沿结果进行分析
% Efficiency=ones(pop_size,11)-P_given/1000;
% figure (1)
% for i=1:pop_size 
% 
%     s(i)=plot(Vin_given,Efficiency(i,:),'LineStyle','--','Marker','o','LineWidth',1.2,'MarkerSize',5);
% 
%     hold on
% Pswitch_off_Matrix(i,1:11)=[Pswitch_off_Matrix_less(i,1:m(i)) Pswitch_off_Matrix_more(i,(m(i)+1):11)];
% P_conduction_Mosfet_Matrix(i,1:11)=[P_conduction_Mosfet_Matrix_less(i,1:m(i)) P_conduction_Mosfet_Matrix_more(i,(m(i)+1):11)];
% P_conduction_Diode_Matrix(i,1:11)=[P_conduction_Diode_Matrix_less(i,1:m(i)) P_conduction_Diode_Matrix_more(i,(m(i)+1):11)];
% P_Tr_Cu_Matrix(i,1:11)=[P_Tr_Cu_Matrix_less(i,1:m(i)) P_Tr_Cu_Matrix_more(i,(m(i)+1):11)];
% P_Lr_Cu_Matrix(i,1:11)=[P_Lr_Cu_Matrix_less(i,1:m(i)) P_Lr_Cu_Matrix_more(i,(m(i)+1):11)];
% 
% P_Tr_Fe_Matrix_less(i,1:m(i))=P_Tr_Fe_t0t2_Matrix_less(i,1:m(i))+P_Tr_Fe_t2t3_Matrix_less(i,1:m(i));
% P_Tr_Fe_Matrix(i,1:11)=[P_Tr_Fe_Matrix_less(i,1:m(i)) P_Tr_Fe_Matrix_more(i,(m(i)+1):11)];
% 
% P_Lr_Fe_Matrix_less(i,1:m(i))=P_Lr_Fe_t0t2_Matrix_less(i,1:m(i))+ P_Lr_Fe_t2t3_Matrix_less(i,1:m(i));
% P_Lr_Fe_Matrix_more(i,(m(i)+1):11)=P_Lr_Fe_t0t2_Matrix_more(i,(m(i)+1):11)+ P_Lr_Fe_t2t3_Matrix_more(i,(m(i)+1):11);  
% P_Lr_Fe_Matrix(i,1:11)=[P_Lr_Fe_Matrix_less(i,1:m(i)) P_Lr_Fe_Matrix_more(i,(m(i)+1):11)];
% P_recovery_Diode_Matrix_more
% 
% end
%     legend(num2str(1),num2str(2),num2str(3),num2str(4),num2str(5),num2str(6),num2str(7),num2str(8),num2str(9),num2str(10),num2str(11));
% 
% %这里的P_total与cal_loss里的平均损耗值应该是一样的，只不过total不参与迭代过程，仅在此用于分析数据
% P_total=Pswitch_off_avg+P_conduction_Mosfet_avg+P_conduction_Diode_avg+P_Tr_Cu_avg+P_Lr_Cu_avg+P_Tr_Fe_avg+P_Lr_Fe_avg+P_recovery_Diode_Matrix_avg;
% Pswitch_off_proportion=Pswitch_off_avg./P_total;
% P_conduction_Mosfet_proportion=P_conduction_Mosfet_avg./P_total;
% P_conduction_Diode_proportion=P_conduction_Diode_avg./P_total;
% P_Tr_Cu_proportion=P_Tr_Cu_avg./P_total;
% P_Lr_Cu_proportion=P_Lr_Cu_avg./P_total;
% P_Tr_Fe_proportion=P_Tr_Fe_avg./P_total;
% P_Lr_Fe_proportion=P_Lr_Fe_avg./P_total;
% P_recovery_Diode_proportion=P_recovery_Diode_Matrix_avg./P_total;
% 
% %%%下面损耗为每个优化个体每个输入电压下的占比
% P_total_vin=Pswitch_off_Matrix+P_conduction_Mosfet_Matrix+P_conduction_Diode_Matrix+...
%             P_Tr_Cu_Matrix+P_Lr_Cu_Matrix+P_Tr_Fe_Matrix+P_Lr_Fe_Matrix+P_recovery_Diode_Matrix_more;
% 
% %开关损耗
% P_switch1=Pswitch_off_Matrix+P_recovery_Diode_Matrix_more;
% P_switch=(Pswitch_off_Matrix+P_recovery_Diode_Matrix_more)./P_total_vin;
% %导通损耗
% P_con1=P_conduction_Mosfet_Matrix+P_conduction_Diode_Matrix;
% P_con=(P_conduction_Mosfet_Matrix+P_conduction_Diode_Matrix)./P_total_vin;
% %磁芯损耗
% P_core1=P_Tr_Fe_Matrix+P_Lr_Fe_Matrix;
% P_core=(P_Tr_Fe_Matrix+P_Lr_Fe_Matrix)./P_total_vin;
% %铜损
% P_copper1=P_Tr_Cu_Matrix+P_Lr_Cu_Matrix;
% P_copper=(P_Tr_Cu_Matrix+P_Lr_Cu_Matrix)./P_total_vin;
% P_histogram=[P_switch P_con P_core P_copper];
% 
% 
% P_total_proportion=[Pswitch_off_proportion P_conduction_Mosfet_proportion P_conduction_Diode_proportion P_Tr_Cu_proportion P_Lr_Cu_proportion P_Tr_Fe_proportion P_Lr_Fe_proportion P_recovery_Diode_proportion]*100;
% P_total_proportion=vpa(P_total_proportion,4)
% vpa(P_total_proportion(:,1),4)%%%关断损耗
% vpa(P_total_proportion(:,2),4)
% vpa(P_total_proportion(:,3),4)
% vpa(P_total_proportion(:,4),4)
% vpa(P_total_proportion(:,5),4)
% vpa(P_total_proportion(:,6),4)
% vpa(P_total_proportion(:,7),4)      
% vpa(P_total_proportion(:,8),4)


 P_given=real(P_given);
%%





end

