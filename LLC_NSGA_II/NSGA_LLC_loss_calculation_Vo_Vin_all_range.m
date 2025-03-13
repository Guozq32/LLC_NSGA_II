function[ Pswitch_off_avg,P_conduction_Mosfet_avg,P_conduction_Diode_avg,P_Tr_Cu_avg,P_Lr_Cu_avg,P_Tr_Fe_avg,P_Lr_Fe_avg, P_recovery_Diode_Matrix_avg,Ploss_Cell ,Pswitch_off_Cell,P_conduction_Mosfet_Cell,P_conduction_Diode_Cell,P_Tr_Cu_Cell,P_Lr_Cu_Cell,Bm_T_t0t2_Cell,P_Tr_Fe_t0t2_Cell,P_Tr_Fe_t2t3_Cell,P_Lr_Fe_t0t2_Cell,P_Lr_Fe_t2t3_Cell,ILs_rms_Cell,Im1_Cell,fs_Cell,Vin_Cell,Vin_Matrix1_Cell ] = NSGA_LLC_loss_calculation_Vo_Vin_all_range(n,Lr,fr,Lm,R,Rac1,Rac2,Rac3,N1,N2,N3,Ae_T,Ae_L,Ve_T,Ve_L,Flag_error)


global Vo Vin_max Vin_min
global Vin_given fs_min fs_max Cr Ln Z0 Z1 wr wn ILs_rms_max Ipeak_max ...
ILs_rms_Matrix_less Ipeak_t0t2_Matrix_less Ipeak_t2t3_Matrix_less Ipeak_Matrix_less fs_Matrix_less  ...
r0_Matrix_less phi0_Matrix_less the0_Matrix_less  ...
r1_Matrix_less phi1_Matrix_less the1_Matrix_less M_Matrix_less Vin_Matrix_less  ...
ILs_rms_Matrix_more Ipeak_Matrix_more fs_Matrix_more  r0_Matrix_more  phi0_Matrix_more the0_Matrix_more  ...
r1_Matrix_more phi1_Matrix_more the1_Matrix_more M_Matrix_more Vin_Matrix_more

global M_max  M_min 

global pop_size tf Ron alpha beta kcof Vo Vin_max   Qrr



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




factor=2;

    for i=1:pop_size
        for j=1:1:10
            if ((n(i)*Vin_given(j))<=Vo)&&((n(i)*Vin_given(j+1))>Vo)
                m(i)=j;
                break

            elseif j==10
                m(i)=11;
            else
                j=j+1;
            end
        end

    end


for i=1:pop_size
    

    if Flag_error(i)==0

            for j=1:m(i)  

                %%%%%%%%%%%%%%switch loss%%%%%%%%%%%%%%
                Im1_Matrix_less(i,j)=Vin_Matrix_less(i,j)*r0_Matrix_less(i,j)*sin(the0_Matrix_less(i,j))/Z0(i);
               
                Pswitch_off_Matrix_less(i,j)=4*Vin_Matrix_less(i,j)*abs(Im1_Matrix_less(i,j))*tf*fs_Matrix_less(i,j)/2; 
                
                %%%%%%%%%%%%%%conduction loss%%%%%%%%%%%%%%
                P_conduction_Mosfet_Matrix_less(i,j)=2*Ron*ILs_rms_Matrix_less(i,j)*ILs_rms_Matrix_less(i,j);
               
                Conduction = @(phi) 0.56 *power(Vin_Matrix_less(i,j)*(r0_Matrix_less(i,j)*sin(phi+the0_Matrix_less(i,j))-r0_Matrix_less(i,j)*sin(the0_Matrix_less(i,j))-M_Matrix_less(i,j)*Lr(i)*phi/Lm(i))/Z0(i)/n(i),0.2233);
                P_conduction_Diode_Matrix_less(i,j)=2*integral(Conduction,0,phi0_Matrix_less(i,j));  
                
                %%%%%%%%%%%%%%  copper loss  %%%%%%%%%%%%%%
                P_Tr_Cu_Matrix_less(i,j)=ILs_rms_Matrix_less(i,j)*ILs_rms_Matrix_less(i,j)*Rac1(i)...
                    +ILs_rms_Matrix_less(i,j)*ILs_rms_Matrix_less(i,j)*Rac2(i)/n(i)/n(i);  
                P_Lr_Cu_Matrix_less(i,j)=ILs_rms_Matrix_less(i,j)*ILs_rms_Matrix_less(i,j)*Rac3(i);                       
               
                %%%%%%%%%%%%%%  magnetic loss  %%%%%%%%%%%%%%
                %%t0t2
                ir_t0_less(i,j)=Vin_Matrix_less(i,j)*r0_Matrix_less(i,j)*sin(the0_Matrix_less(i,j))/Z0(i);
                ir_t2_less(i,j)=Vin_Matrix_less(i,j)*r0_Matrix_less(i,j)*sin(the0_Matrix_less(i,j)+phi0_Matrix_less(i,j))/Z0(i);
                deta_B_T_t0t2_less(i,j)=Lm(i)*(ir_t2_less(i,j)-ir_t0_less(i,j))/N1(i)/Ae_T(i);
             
                dB_T_t0t2_less(i,j)= Vo/N2(i)/Ae_T(i);
                P_Tr_Fe_t0t2_Matrix_less(i,j)=(2*fs_Matrix_less(i,j)*phi0_Matrix_less(i,j)/wr(i))*kcof*Ve_T(i)... 
                    *power(dB_T_t0t2_less(i,j),alpha)*power( deta_B_T_t0t2_less(i,j),(beta-alpha));
            
                %%t2t3
               
                ir_t3_less(i,j)=Vin_Matrix_less(i,j)*r1_Matrix_less(i,j)*sin(the1_Matrix_less(i,j)+phi1_Matrix_less(i,j))/Z1(i);
            
                deta_i_t2t3_less(i,j)=max(abs(Ipeak_t2t3_Matrix_less(i,j)-ir_t2_less(i,j)),abs(Ipeak_t2t3_Matrix_less(i,j)-ir_t3_less(i,j)));
                deta_B_T_t2t3_less(i,j)=Lm(i)*(deta_i_t2t3_less(i,j))/N1(i)/Ae_T(i);
                
                fun1_less = @(phi) power(abs(Lm(i)*Vin_Matrix_less(i,j)*r1_Matrix_less(i,j)*wn(i)*cos(phi+the1_Matrix_less(i,j))/Z1(i)/N1(i)/Ae_T(i)),alpha);
                
                P_Tr_Fe_t2t3_Matrix_less(i,j)=2*fs_Matrix_less(i,j)...
                *power(deta_B_T_t2t3_less(i,j),beta-alpha)*integral(fun1_less,0,phi1_Matrix_less(i,j))*kcof*Ve_T(i)/wn(i);
            
            
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
            end        

            if m(i)~= 11
            
       

            for j=(m(i)+1):11  
                %%%%%%%%%%%%%%switch loss%%%%%%%%%%%%%%       

                Im1_Matrix_more(i,j)=Vin_Matrix_more(i,j)*r0_Matrix_more(i,j)*sin(the0_Matrix_more(i,j)+phi0_Matrix_more(i,j))/Z0(i);
               
                Pswitch_off_Matrix_more(i,j)=4*Vin_Matrix_more(i,j)*abs(Im1_Matrix_more(i,j))*tf*fs_Matrix_more(i,j)/2;                 
                
                %%%%%%%%%%%%%%conduction loss%%%%%%%%%%%%%%

                P_conduction_Mosfet_Matrix_more(i,j)=2*Ron*ILs_rms_Matrix_more(i,j)*ILs_rms_Matrix_more(i,j);
                
                Conduction_more_t0t2 = @(phi) 0.56*power(Vin_Matrix_more(i,j)*(r0_Matrix_more(i,j)*sin(phi+the0_Matrix_more(i,j))-r0_Matrix_more(i,j)*sin(the0_Matrix_more(i,j))-M_Matrix_more(i,j)*Lr(i)*phi/Lm(i))/Z0(i)/n(i),0.2233);
                Conduction_more_t2t3 = @(phi) 0.56*power(Vin_Matrix_more(i,j)*(r1_Matrix_more(i,j)*sin(phi+the1_Matrix_more(i,j))-r0_Matrix_more(i,j)*sin(the0_Matrix_more(i,j))-M_Matrix_more(i,j)*Lr(i)*(phi+phi0_Matrix_more(i,j))/Lm(i))/Z0(i)/n(i),0.2233);
                P_conduction_Diode_Matrix_more(i,j)=2*integral(Conduction_more_t0t2,0,phi0_Matrix_more(i,j))+2*integral(Conduction_more_t2t3,0,phi1_Matrix_more(i,j));   
             
               
                %%%%%%%%%%%%%%Diode reverse recovery loss%%%%%%%%%%%%%%
                
                P_recovery_Diode_Matrix_more(i,j)=4*Vo*Qrr*fs_Matrix_more(i,j);     

                %%%%%%%%%%%%%%  copper loss  %%%%%%%%%%%%%%
                P_Tr_Cu_Matrix_more(i,j)=ILs_rms_Matrix_more(i,j)*ILs_rms_Matrix_more(i,j)*Rac1(i)+ILs_rms_Matrix_more(i,j)*ILs_rms_Matrix_more(i,j)*Rac2(i)/n(i)/n(i);   
                
                P_Lr_Cu_Matrix_more(i,j)=ILs_rms_Matrix_more(i,j)*ILs_rms_Matrix_more(i,j)*Rac3(i);                     
               
                %%%%%%%%%%%%%%  ´ÅËð  %%%%%%%%%%%%%%

             
                deta_B_T_more(i,j)= Vo/N2(i)/2/fs_Matrix_more(i,j)/Ae_T(i);       
                dB_T_more(i,j)= Vo/N2(i)/Ae_T(i);
                P_Tr_Fe_Matrix_more(i,j)=(2*fs_Matrix_more(i,j)*pi/wr(i))*((power(dB_T_more(i,j),alpha)*power(deta_B_T_more(i,j),beta-alpha))*kcof*Ve_T(i));             
            
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

            end
            else
                m(i)=11;
            end

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


            %%%%%%%%%%%%%-----Average loss calculation----%%%%%%%%%%%%%%
            
            Pswitch_off_avg(i)=(sum(Pswitch_off_Matrix_less(i,1:m(i)))+sum(Pswitch_off_Matrix_more(i,(m(i)+1):11)))/11;  
            P_conduction_Mosfet_avg(i)=(sum(P_conduction_Mosfet_Matrix_less(i,1:m(i)))+sum(P_conduction_Mosfet_Matrix_more(i,(m(i)+1):11)))/11 ;
            P_conduction_Diode_avg(i)=(sum(P_conduction_Diode_Matrix_less(i,1:m(i)))+sum(P_conduction_Diode_Matrix_more(i,(m(i)+1):11)))/11;
            
            P_recovery_Diode_Matrix_avg(i)=sum(P_recovery_Diode_Matrix_more(i,(m(i)+1):11))/11;
            P_Tr_Cu_avg(i)=(sum(P_Tr_Cu_Matrix_less(i,1:m(i)))+sum(P_Tr_Cu_Matrix_more(i,(m(i)+1):11)))/11;  
            P_Lr_Cu_avg(i)=(sum(P_Lr_Cu_Matrix_less(i,1:m(i)))+sum(P_Lr_Cu_Matrix_more(i,(m(i)+1):11)))/11;  
          
            P_Tr_Fe_avg(i)=(sum(P_Tr_Fe_t0t2_Matrix_less(i,1:m(i)))+sum(P_Tr_Fe_t2t3_Matrix_less(i,1:m(i)))+sum(P_Tr_Fe_Matrix_more(i,(m(i)+1):11)))/11;
      
           P_Lr_Fe_t0t2_avg(i)=(sum(P_Lr_Fe_t0t2_Matrix_less(i,1:m(i)))+sum(P_Lr_Fe_t0t2_Matrix_more(i,(m(i)+1):11)))/11;            
           P_Lr_Fe_t2t3_avg(i)=(sum(P_Lr_Fe_t2t3_Matrix_less(i,1:m(i)))+sum(P_Lr_Fe_t2t3_Matrix_more(i,(m(i)+1):11)))/11;
                  
            P_Lr_Fe_avg(i)=P_Lr_Fe_t0t2_avg(i)+P_Lr_Fe_t2t3_avg(i);
    end
end



%%%%%%%%%%%%%-----The conversion parameters are stored as cell-----%%%%%%%%%%%%%%
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



end

