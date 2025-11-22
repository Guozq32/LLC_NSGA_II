function [k3_cell,Ipeak_cell,ILs_rms_cell,Irms_p_cell,Irms_s_cell,Pswitch_on_cell,Pswitch_off_cell,Pswitch_on_diode_cell,P_Tr_Cu_cell,P_Tr_Fe_cell,P_Lr_Cu_cell,P_Lr_Fe_cell,Ploss_cell,Ploss_percentage_cell,Efficiency_cell] = transform(k3,Ipeak,ILs_rms,Irms_p,Irms_s,Pswitch_on,Pswitch_off,Pswitch_on_diode,P_Tr_Cu,P_Tr_Fe,P_Lr_Cu,P_Lr_Fe,Ploss,Ploss_percentage,Efficiency)
global pop_size

for i=1:pop_size 
    k3_cell(i)={k3(i,:)};
    Ipeak_cell(i)={Ipeak(i,:)};
    ILs_rms_cell(i)={ILs_rms(i,:)};
    Irms_p_cell(i)={Irms_p(i,:)};
    Irms_s_cell(i)={Irms_s(i,:)};
    Pswitch_on_cell(i)={Pswitch_on(i,:)};
    Pswitch_off_cell(i)={Pswitch_off(i,:)};
    Pswitch_on_diode_cell(i)={Pswitch_on_diode(i,:)};
    P_Tr_Cu_cell(i)={P_Tr_Cu(i,:)};
    P_Tr_Fe_cell(i)={P_Tr_Fe(i,:)};
    P_Lr_Cu_cell(i)={P_Lr_Cu(i,:)};
    P_Lr_Fe_cell(i)={P_Lr_Fe(i,:)};
    Ploss_cell(i)={Ploss(i,:)};
    Ploss_percentage_cell(i)={Ploss_percentage(i,:)};
    Efficiency_cell(i)={Efficiency(i,:)};
end
k3_cell=k3_cell';
Ipeak_cell=Ipeak_cell';
ILs_rms_cell=ILs_rms_cell';
Irms_p_cell=Irms_p_cell';
Irms_s_cell=Irms_s_cell';
Pswitch_on_cell=Pswitch_on_cell';
Pswitch_off_cell=Pswitch_off_cell';
Pswitch_on_diode_cell= Pswitch_on_diode_cell';
P_Tr_Cu_cell=P_Tr_Cu_cell';
P_Tr_Fe_cell=P_Tr_Fe_cell';
P_Lr_Cu_cell=P_Lr_Cu_cell';
P_Lr_Fe_cell=P_Lr_Fe_cell';
Ploss_cell=Ploss_cell';
Ploss_percentage_cell=Ploss_percentage_cell';
Efficiency_cell=Efficiency_cell';
end