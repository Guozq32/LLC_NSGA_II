function [x,Para,Para_cell,ff,err] = Cal_Loss(x) 
global Po pop_size Vo Vin_min Vin_max rou miu ds diso xl xu  Ae Aw Ve AeAw h hw Icm Sm DL tdead Coss eta tf alpha beta kcof Ron%%rou表示铜的电阻率，miu表示真空磁导率，表示DL对应利兹线线径，ds对应利兹线单股线直径, Sm1表示原边利兹线股数，Sm2表示副边利兹线股数

global fs_min fs_max Cr Ln Z0 Z1 wr wn ILs_rms_max Ipeak_max ...
ILs_rms_Matrix_less Ipeak_t0t2_Matrix_less Ipeak_t2t3_Matrix_less Ipeak_Matrix_less fs_Matrix_less  ...
r0_Matrix_less phi0_Matrix_less the0_Matrix_less  ...
r1_Matrix_less phi1_Matrix_less the1_Matrix_less M_Matrix_less Vin_Matrix_less  ...
ILs_rms_Matrix_more Ipeak_Matrix_more fs_Matrix_more  r0_Matrix_more  phi0_Matrix_more the0_Matrix_more  ...
r1_Matrix_more phi1_Matrix_more the1_Matrix_more M_Matrix_more Vin_Matrix_more

global M_max  M_min 




n=x(:,1);
Lr=x(:,2);
fr=x(:,3);
K=x(:,4);
Lm=x(:,2).*x(:,4);



M_min  = zeros(pop_size,1);
M_max  = zeros(pop_size,1);

fs_min  = zeros(pop_size,1);
fs_max  = zeros(pop_size,1);
Cr      = zeros(pop_size,1);
Ln      = zeros(pop_size,1);
Z0      = zeros(pop_size,1);
Z1      = zeros(pop_size,1);
wr      = zeros(pop_size,1); 
wn      = zeros(pop_size,1);
ILs_rms_max = zeros(pop_size,1);
Ipeak_max   = zeros(pop_size,1);


J=500e4; 
Kw_Tr=0.13;
Kw_Lr=0.13;

for i=1:pop_size
    if n(i)<xl(1)||isnan(n(i))==1
        n(i)=xl(1);
    elseif n(i)>xu(1)
        n(i)=xu(1);
    else
    end
end
n_1=n;
for i=1:pop_size
    if Lr(i)<xl(2)||isnan(Lr(i))==1
        Lr(i)=xl(2);
    elseif Lr(i)>xu(2)
        Lr(i)=xu(2);
    else
    end
end
Lr_1=Lr;
for i=1:pop_size
    if fr(i)<xl(3)||isnan(fr(i))==1
        fr(i)=xl(3);
    elseif fr(i)>xu(3)
        fr(i)=xu(3);
    else
    end
end

Bm_T=x(:,5);
Bm_L=x(:,6);
for i=1:pop_size
    if Bm_T(i)<xl(5)||isnan(Bm_T(i))==1
        Bm_T(i)=xl(5);
    elseif Bm_T(i)>xu(5)
        Bm_T(i)=xu(5);
    else
    end
end
for i=1:pop_size
    if Bm_L(i)<xl(6)||isnan(Bm_L(i))==1
        Bm_L(i)=xl(6);
    elseif Bm_L(i)>xu(6)
        Bm_L(i)=xu(6);
    else
    end
end

R=Vo*Vo/Po;
[Flag_error,error_1,Flag_PN,Flag_OPO,Flag_PON,Flag_NOP,Flag_ZVS,Flag_M_less,Flag_M_more ]= NSGA_LLC_steady_cal_fs_all_range(n,Lr,fr,Lm,R);



Tr_index=zeros(pop_size,1);
h_T=zeros(pop_size,1);
hw_T=zeros(pop_size,1);
Ve_T=zeros(pop_size,1);
Ae_T=zeros(pop_size,1);
Aw_T=zeros(pop_size,1);
AeAw_Tr=zeros(pop_size,1);
Sm1=zeros(pop_size,1);
DL1=zeros(pop_size,1);
Sm2=zeros(pop_size,1);
DL2=zeros(pop_size,1);
Irms_s_max=zeros(pop_size,1);     
Nl1=zeros(pop_size,1);
Nl2=zeros(pop_size,1);
diso_T=diso*ones(pop_size,1);
N2=zeros(pop_size,1);
N1=zeros(pop_size,1);
N2_floor=zeros(pop_size,1);
N1_floor=zeros(pop_size,1);
N2_ceil=zeros(pop_size,1);
N1_ceil=zeros(pop_size,1);
N_Matrix=zeros(pop_size,4);
N_min=zeros(pop_size,1);
n_Matrix=zeros(pop_size,4);
N1_Martix=zeros(pop_size,4);
N2_Martix=zeros(pop_size,4);
m1=zeros(pop_size,1);
m2=zeros(pop_size,1);
Lr_index=zeros(pop_size,1);
h_L=zeros(pop_size,1);
hw_L=zeros(pop_size,1);
Ve_L=zeros(pop_size,1);
Ae_L=zeros(pop_size,1);
Aw_L=zeros(pop_size,1);
AeAw_L=zeros(pop_size,1);
DL3=zeros(pop_size,1);
Sm3=zeros(pop_size,1);
N3=zeros(pop_size,1);
Nl3=zeros(pop_size,1);
m3=zeros(pop_size,1);


for i=1:pop_size
    
    if Flag_error(i)==0  

          AeAw_Tr(i)=(1+eta)*Po/(4*eta*Bm_T(i)*J*Kw_Tr*fs_min(i));

            if AeAw_Tr(i)<=AeAw(1)
                Tr_index(i)=1;
                AeAw_Tr(i)=AeAw(1);
            elseif AeAw_Tr(i)>AeAw(1) && AeAw_Tr(i)<=AeAw(2)
                Tr_index(i)=2;
                AeAw_Tr(i)=AeAw(2);
            elseif AeAw_Tr(i)>AeAw(2) && AeAw_Tr(i)<=AeAw(3)
                Tr_index(i)=3;
                AeAw_Tr(i)=AeAw(3);
            elseif AeAw_Tr(i)>AeAw(3) && AeAw_Tr(i)<=AeAw(4)
                Tr_index(i)=4;
                AeAw_Tr(i)=AeAw(4);
            elseif AeAw_Tr(i)>AeAw(4) && AeAw_Tr(i)<=AeAw(5)
                Tr_index(i)=5;
                AeAw_Tr(i)=AeAw(5);
            else
                Tr_index(i)=6;
                AeAw_Tr(i)=AeAw(6);
            end
            h_T(i)=h(Tr_index(i));
            hw_T(i)=hw(Tr_index(i));
            Ve_T(i)=Ve(Tr_index(i));
            Ae_T(i)=Ae(Tr_index(i));
            Aw_T(i)=Aw(Tr_index(i));

            if ILs_rms_max(i)<=Icm(1)
                Sm1(i)=Sm(1);
                DL1(i)=DL(1);
            elseif ILs_rms_max(i)>Icm(1) && ILs_rms_max(i)<=Icm(2)
                Sm1(i)=Sm(2);
                DL1(i)=DL(2);
            elseif ILs_rms_max(i)>Icm(2) && ILs_rms_max(i)<=Icm(3)
                Sm1(i)=Sm(3);
                DL1(i)=DL(3);
            elseif ILs_rms_max(i)>Icm(3) && ILs_rms_max(i)<=Icm(4)
                Sm1(i)=Sm(4);
                DL1(i)=DL(4);
            elseif ILs_rms_max(i)>Icm(4) && ILs_rms_max(i)<=Icm(5)
                Sm1(i)=Sm(5);
                DL1(i)=DL(5);
            elseif ILs_rms_max(i)>Icm(5) && ILs_rms_max(i)<=Icm(6)
                Sm1(i)=Sm(6);
                DL1(i)=DL(6);
            elseif ILs_rms_max(i)>Icm(6) && ILs_rms_max(i)<=Icm(7)
                Sm1(i)=Sm(7);
                DL1(i)=DL(7);
            elseif ILs_rms_max(i)>Icm(7) && ILs_rms_max(i)<=Icm(8)
                Sm1(i)=Sm(8);
                DL1(i)=DL(8);
            elseif ILs_rms_max(i)>Icm(8) && ILs_rms_max(i)<=Icm(9)
                Sm1(i)=Sm(9);
                DL1(i)=DL(9);
            elseif ILs_rms_max(i)>Icm(9) && ILs_rms_max(i)<=Icm(10)
                Sm1(i)=Sm(10);
                DL1(i)=DL(10);
            elseif ILs_rms_max(i)>Icm(10) && ILs_rms_max(i)<=Icm(11)
                Sm1(i)=Sm(11);
                DL1(i)=DL(11);
            elseif ILs_rms_max(i)>Icm(11) && ILs_rms_max(i)<=Icm(12)
                Sm1(i)=Sm(12);
                DL1(i)=DL(12);
            else
                Sm1(i)=Sm(13);
                DL1(i)=DL(13);
            end


        Irms_s_max(i)=ILs_rms_max(i)/n(i);

            if Irms_s_max(i)<=Icm(1)
                Sm2(i)=Sm(1);
                DL2(i)=DL(1);
            elseif Irms_s_max(i)>Icm(1) && Irms_s_max(i)<=Icm(2)
                Sm2(i)=Sm(2);
                DL2(i)=DL(2);
            elseif Irms_s_max(i)>Icm(2) && Irms_s_max(i)<=Icm(3)
                Sm2(i)=Sm(3);
                DL2(i)=DL(3);
            elseif Irms_s_max(i)>Icm(3) && Irms_s_max(i)<=Icm(4)
                Sm2(i)=Sm(4);
                DL2(i)=DL(4);
            elseif Irms_s_max(i)>Icm(4) && Irms_s_max(i)<=Icm(5)
                Sm2(i)=Sm(5);
                DL2(i)=DL(5);
            elseif Irms_s_max(i)>Icm(5) && Irms_s_max(i)<=Icm(6)
                Sm2(i)=Sm(6);
                DL2(i)=DL(6);
            elseif Irms_s_max(i)>Icm(6) && Irms_s_max(i)<=Icm(7)
                Sm2(i)=Sm(7);
                DL2(i)=DL(7);
            elseif Irms_s_max(i)>Icm(7) && Irms_s_max(i)<=Icm(8)
                Sm2(i)=Sm(8);
                DL2(i)=DL(8);
            elseif Irms_s_max(i)>Icm(8) && Irms_s_max(i)<=Icm(9)
                Sm2(i)=Sm(9);
                DL2(i)=DL(9);
            elseif Irms_s_max(i)>Icm(9) && Irms_s_max(i)<=Icm(10)
                Sm2(i)=Sm(10);
                DL2(i)=DL(10);
            elseif Irms_s_max(i)>Icm(10) && Irms_s_max(i)<=Icm(11)
                Sm2(i)=Sm(11);
                DL2(i)=DL(11);
            elseif Irms_s_max(i)>Icm(11) && Irms_s_max(i)<=Icm(12)
                Sm2(i)=Sm(12);
                DL2(i)=DL(12);
            else
                Sm2(i)=Sm(13);
                DL2(i)=DL(13);
            end

        Nl1(i)=floor(h_T(i)/DL1(i));  
        Nl2(i)=floor(h_T(i)/DL2(i));
        N2(i)=Vo/(4*Bm_T(i)*Ae_T(i)*fs_min(i)); 
        N1(i)=N2(i)/n(i);
        N2_floor(i)=floor(N2(i));
        N1_floor(i)=floor(N1(i));
        N2_ceil(i)=ceil(N2(i));
        N1_ceil(i)=ceil(N1(i));
        N_Matrix(i,1:4)=[abs(n(i)-(N2_floor(i)/N1_floor(i))),abs(n(i)-(N2_floor(i)/N1_ceil(i))),abs(n(i)-(N2_ceil(i)/N1_floor(i))),abs(n(i)-(N2_ceil(i)/N1_ceil(i)))];
        [~,N_min(i)]=min(N_Matrix(i,1:4));
        n_Matrix(i,1:4)=[N2_floor(i)/N1_floor(i),N2_floor(i)/N1_ceil(i),N2_ceil(i)/N1_floor(i),N2_ceil(i)/N1_ceil(i)];
        N1_Martix(i,1:4)=[N1_floor(i),N1_ceil(i),N1_floor(i),N1_ceil(i)];
        N2_Martix(i,1:4)=[N2_floor(i),N2_floor(i),N2_ceil(i),N2_ceil(i)];
        N1(i)=N1_Martix(i,N_min(i));
        N2(i)=N2_Martix(i,N_min(i));
        m1(i)=ceil(N1(i)/Nl1(i));
        m2(i)=ceil(N2(i)/Nl2(i));
 
        AeAw_L(i)=Lr(i)*Ipeak_max(i)*ILs_rms_max(i)/(Kw_Lr*J*Bm_L(i));
            if AeAw_L(i)<=AeAw(1)
                Lr_index(i)=1;
                AeAw_L(i)=AeAw(1);
            elseif AeAw_L(i)>AeAw(1) && AeAw_L(i)<=AeAw(2)
                Lr_index(i)=2;
                AeAw_L(i)=AeAw(2);
            elseif AeAw_L(i)>AeAw(2) && AeAw_L(i)<=AeAw(3)
                Lr_index(i)=3;
                AeAw_L(i)=AeAw(3);
            elseif AeAw_L(i)>AeAw(3) && AeAw_L(i)<=AeAw(4)
                Lr_index(i)=4;
                AeAw_L(i)=AeAw(4);
            elseif AeAw_L(i)>AeAw(4) && AeAw_L(i)<=AeAw(5)
                Lr_index(i)=5;
                AeAw_L(i)=AeAw(5);
            else
                Lr_index(i)=6;
                AeAw_L(i)=AeAw(6);
            end
            h_L(i)=h(Lr_index(i));
            hw_L(i)=hw(Lr_index(i));
            Ve_L(i)=Ve(Lr_index(i));
            Ae_L(i)=Ae(Lr_index(i));
            Aw_L(i)=Aw(Lr_index(i));
            DL3(i)=DL1(i); 
            Sm3(i)=Sm1(i);
            N3(i)=ceil(Lr(i)*Ipeak_max(i)/(Bm_L(i)*Ae_L(i)));
            Nl3(i)=floor(h_L(i)/DL3(i));
            m3(i)=ceil(N3(i)/Nl3(i));
    end
end

Km_T=zeros(pop_size,1);
Km_L=zeros(pop_size,1);
lp=zeros(pop_size,1);      
lw=zeros(pop_size,1);
ls=zeros(pop_size,1); 
lr=zeros(pop_size,1);
Rdc1=zeros(pop_size,1);
Rdc2=zeros(pop_size,1); 
Rdc3=zeros(pop_size,1);
Ddelta=zeros(pop_size,1);
kxi1=zeros(pop_size,1);
kxi2=zeros(pop_size,1); 
kxi3=zeros(pop_size,1);
Rac1=zeros(pop_size,1);
Rac2=zeros(pop_size,1); 
Rac3=zeros(pop_size,1);
Fac1=zeros(pop_size,1);
Fac2=zeros(pop_size,1);
Fac3=zeros(pop_size,1);
error_3=zeros(pop_size,1);
Flag_error_3=zeros(pop_size,1);
factor=1.5;        
for i=1:pop_size

    Km_T(i)=(pi*(ds/2)^2*Sm1(i)*N1(i)+pi*(ds/2)^2*Sm2(i)*N2(i))/Aw_T(i); 
    Km_L(i)=(pi*(ds/2)^2*Sm3(i)*N3(i))/Aw_L(i);
    if Km_T(i)>0.205 || Km_L(i)>0.205
        error_3(i)=100;
        Flag_error_3(i)=1;
    end
            lp(i)=sqrt(Ae_T(i)/pi)+(m1(i)*DL1(i)+diso*(m1(i)-1))*0.5; 
            lw(i)=sqrt(Ae_T(i)/pi)+m1(i)*DL1(i)+diso*m1(i); 
            ls(i)=lw(i)+(m2(i)*DL2(i)+diso.*(m2(i)-1))*0.5; 
            lr(i)=sqrt(Ae_L(i)/pi)+(m3(i)*DL3(i)+diso*(m3(i)-1))*0.5; 
            Rdc1(i)=rou*2*pi*lp(i)*N1(i)/(pi*(ds/2)^2*Sm1(i)); 
            Rdc2(i)=rou*2*pi*ls(i)*N2(i)/(pi*(ds/2)^2*Sm2(i)); 
            Rdc3(i)=rou*2*pi*lr(i)*N3(i)/(pi*(ds/2)^2*Sm3(i));  
             Ddelta(i)=2^0.5/sqrt(2*pi*fr(i)*miu/rou); 
    
            kxi1(i)=0.5*ds/Ddelta(i)*sqrt(pi*Nl1(i)*sqrt(Sm1(i)*pi)*ds*0.5/h_T(i)); 
            kxi2(i)=0.5*ds/Ddelta(i)*sqrt(pi*Nl2(i)*sqrt(Sm2(i)*pi)*ds*0.5/h_T(i));  
            kxi3(i)=0.5*ds/Ddelta(i)*sqrt(pi*Nl3(i)*sqrt(Sm3(i)*pi)*ds*0.5/h_L(i)); 
    
            k11(i)=vpa(sinh(vpa(2*kxi1(i),32))+sin(vpa(2*kxi1(i),32)),32);
            k12(i)=vpa(cosh(vpa(2*kxi1(i),32))-cos(vpa(2*kxi1(i),32)),32);
            k13(i)=vpa(sinh(vpa(kxi1(i),32))-sin(vpa(kxi1(i),32)),32);
            k14(i)=vpa(cosh(vpa(kxi1(i),32))+cos(vpa(kxi1(i),32)),32);
        
            Fac1(i)=vpa(kxi1(i))*(vpa(k11(i)/k12(i))+2/3*(m1(i)^2-1)*vpa(k13(i)/k14(i)));
            Fac2(i)=kxi2(i)*((sinh(2*kxi2(i))+sin(2*kxi2(i)))/(cosh(2*kxi2(i))-cos(2*kxi2(i)))+2/3*(m2(i)^2-1)*(sinh(kxi2(i))-sin(kxi2(i)))/(cosh(kxi2(i))+cos(kxi2(i))));
            Fac3(i)=kxi3(i)*((sinh(2*kxi3(i))+sin(2*kxi3(i)))/(cosh(2*kxi3(i))-cos(2*kxi3(i)))+2/3*(m3(i)^2-1)*(sinh(kxi3(i))-sin(kxi3(i)))/(cosh(kxi3(i))+cos(kxi3(i))));
   
            Rac1(i)=Rdc1(i)*Fac1(i)*factor;  
            Rac2(i)=Rdc2(i)*Fac2(i)*factor;  
            Rac3(i)=Rdc3(i)*Fac3(i)*factor;  
end



[ Pswitch_off_avg,P_conduction_Mosfet_avg,P_conduction_Diode_avg,P_Tr_Cu_avg,P_Lr_Cu_avg,P_Tr_Fe_avg,P_Lr_Fe_avg, P_recovery_Diode_Matrix_avg,Ploss_Cell ,Pswitch_off_Cell,P_conduction_Mosfet_Cell,P_conduction_Diode_Cell,P_Tr_Cu_Cell,P_Lr_Cu_Cell,Bm_T_t0t2_Cell,P_Tr_Fe_t0t2_Cell,P_Tr_Fe_t2t3_Cell,P_Lr_Fe_t0t2_Cell,P_Lr_Fe_t2t3_Cell,ILs_rms_Cell,Im1_Cell,fs_Cell,Vin_Cell,Vin_Matrix1_Cell ] = NSGA_LLC_loss_calculation_Vo_Vin_all_range(n,Lr,fr,Lm,R,Rac1,Rac2,Rac3,N1,N2,N3,Ae_T,Ae_L,Ve_T,Ve_L,Flag_error);
1




error=max(error_1,error_3);

Ploss_avg=Pswitch_off_avg+P_conduction_Mosfet_avg+P_conduction_Diode_avg+P_Tr_Cu_avg+P_Lr_Cu_avg+P_Tr_Fe_avg+P_Lr_Fe_avg+P_recovery_Diode_Matrix_avg;


ff=[Ploss_avg Ve_T+Ve_L];
x=[n,Lr,fr,K,Bm_T,Bm_L];
Para=[N1 N2 N3 Nl1 Nl2 Nl3 m1 m2 m3 Rac1 ...
    Rac2 Rac3 Sm1 Sm2 Sm3 Km_T Km_L Ve_T Ve_L Tr_index ...
    Lr_index Lm  Cr Ploss_avg Pswitch_off_avg ...
    P_conduction_Mosfet_avg P_conduction_Diode_avg P_Tr_Cu_avg P_Lr_Cu_avg P_Tr_Fe_avg ...
    P_Lr_Fe_avg  P_recovery_Diode_Matrix_avg Ipeak_max ILs_rms_max  fs_max fs_min  error];

Para_cell=[Ploss_Cell Pswitch_off_Cell P_conduction_Mosfet_Cell P_conduction_Diode_Cell P_Tr_Cu_Cell P_Lr_Cu_Cell P_Tr_Fe_t0t2_Cell P_Tr_Fe_t2t3_Cell P_Lr_Fe_t0t2_Cell P_Lr_Fe_t2t3_Cell ILs_rms_Cell Im1_Cell fs_Cell Vin_Cell Vin_Matrix1_Cell Bm_T_t0t2_Cell];

err=error;


end



