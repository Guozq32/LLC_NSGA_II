clear all
clc
close all
% format long
global Po Vo Vin_min Vin_max rou miu ds diso Rdson alpha beta kcof etac etam pop_size pm V M xl xu Ae Aw Ve AeAw h hw Icm Sm DL tdead Coss eta tf Ron  Qrr



global Vin_given fs_min fs_max Cr Ln Z0 Z1 wr wn ILs_rms_max Ipeak_max ...
ILs_rms_Matrix_less Ipeak_t0t2_Matrix_less Ipeak_t2t3_Matrix_less Ipeak_Matrix_less fs_Matrix_less  ...
r0_Matrix_less phi0_Matrix_less the0_Matrix_less  ...
r1_Matrix_less phi1_Matrix_less the1_Matrix_less M_Matrix_less Vin_Matrix_less  ...
ILs_rms_Matrix_more Ipeak_Matrix_more fs_Matrix_more  r0_Matrix_more  phi0_Matrix_more the0_Matrix_more  ...
r1_Matrix_more phi1_Matrix_more the1_Matrix_more M_Matrix_more Vin_Matrix_more

global M_max M_min

etac = 0.3;                  % distribution index for crossover 
etam = 0.3;                  % distribution index for mutation / mutation constant 
M=2;                         % Number of optimization objectives      
pop_size=200;           % Population size 
gen_max=3;            % MAx number of generations - stopping criteria 

%%%% In this programmer£¬ n is defined as vo/vin
Vo=200;
Vin_min=150;
Vin_max=250;
Po=1000;
rou=2.12e-8;	%%Electrical resistivity of copper  0.0185¦¸¡¤mm2/m->@25¡æ;  0.0212¦¸¡¤mm2/m->@70¡æ;
miu=4*pi*1e-7;  %%Permeability of free space
Qgd_p=7e-9;                        %gate leakage charge                    
Qgs_p=8e-9;                        %gate source charge                    
Ciss_p=640e-12;                    %gate source charge                    
Crss_p=2.3e-12;                      %Inverse conductive capacitance                  
Vsp_p=6;                            %Miller platform voltage      
VGS_th_p=2.8;                       %gate threshold voltage                 
RG_int_p=6;                       %Internal gate resistance               
Vcc=15;                             %High driving voltage level        
Vee=-5;                             %Low driving voltage level          
Rdriver_off=2;                    %driving resistance    
RG_off_p=Rdriver_off+RG_int_p;      %Turn-off resistance
Cgs_p=Ciss_p-Crss_p;                %Internal gate-source capacitance
Qth_off_p=(VGS_th_p-Vee)*Cgs_p;
QG_off_p=Qgd_p+Qgs_p-Qth_off_p;
IG_off_p=(Vsp_p-Vee)/RG_off_p;
tf=QG_off_p/IG_off_p;           %turn-off time

tdead=60e-9;  %%% Dead time
Qrr= 20e-9;   %%% Diode reverse recovery charge,refer to the manual for details
Coss=260e-12; %% MOSFET parasitic capacitance
eta=1;    
ds=0.1e-3;  %%Diameter of each strand in Litz wire, unit in meters
diso=0.1e-3; %%Insulation distance between each winding layer
Ron=0.08;% MOSFET on-resistance
alpha=1.658963082164934; 
ki=0.123934484240003;
beta=2.410283969308211;
fun=@(theta) power(abs(cos(theta)),alpha)*power(2,beta-alpha);
q=integral(fun,0,2*pi);
kcof=ki/(power(2*pi,alpha-1)*q);


%     1      2       3       4       5       6
%  PQ26/25  PQ32/20 PQ32/30 PQ35/35 PQ40/40 PQ50/50
Ae = [120	169     167     190     202     328]*1e-6;    
Aw = [84.53	80.50	149.1	220     325.98	433.20]*1e-6; 
AeAw=[10143 13604.5 24899.7 41800   65846.95 142089.6]*1e-12; 
Ve = [6530  9440    12500   16300   20500   37100]*1e-9;  
h  = 0.7*[16.1  11.5    21.3    25      29.5    36.1]*1e-3; 
hw = 0.7*[5.3   7       7       8.8     11.1    12]*1e-3;   
%     1     2   3    4   5    6   7     8    9    10    11  12  13
Icm= [3.95 4.3  4.7  5.1 5.5  5.9 6.3   7.05 7.85 8.65  9.8	11	11.8]; 
Sm = [100  110  120  130 140  150 160   180  200  220   250	280 300]; 
DL = [1.4  1.47 1.53 1.6 1.66 1.71 1.77 1.88 1.98 2.08  2.21 2.34 2.42]*1e-3;

%   1   2      3     4      5     6       %   1 2   3  4    5    6    
xl=[0.82 8e-6  60e3  2  0.13 0.13];   %%x=[n,Lr,fr,K,  Bm_T,Bm_L]  
xu=[1.2 50e-6 100e3  7.5  0.25 0.25];    %%x=[n,Lr,fr,K,  Bm_T,Bm_L]
V=size(xl,2)
pm=1/V;                     % Mutation Probability 



xl_temp=repmat(xl, pop_size,1);  %% xl _ temp is pop _ size ¡Á V dimension. The number of rows of the matrix represents the number of populations, and the number of columns represents the number of optimization objective variables.
xu_temp=repmat(xu, pop_size,1);
x = xl_temp+((xu_temp-xl_temp).*rand(pop_size,V)); 
[x] = Re_random(x,xl,xu) ;


%% Evaluate objective function
[x,Para,Para_cell,ff,err] = Cal_Loss(x);%ff two column optimization target, err violation value


error_norm=normalisation(err);                     % Normalisation of the constraint violation  
population_init=[x ff error_norm]; 
[population,initial_data_Para_sorted,initial_data_Para_cell_sorted,front]=NDS_CD_cons(population_init,Para,Para_cell);    % Non domination Sorting on initial population 
population_initial=population;
population_front1=population(find(population(:,V+M+2)==1),:);
population_1_feasible=population(find(population(:,V+M+1)==0),:);
population_1_infeasible=population(find(population(:,V+M+1)~=0),:);
parent_data_Para=initial_data_Para_sorted;
parent_data_Para_cell=initial_data_Para_cell_sorted;
%% Generation Starts
for gen_count=1:gen_max
    gen_count 
parent_selected=tour_selection(population);                     
child_offspring  = genetic_operator(parent_selected(:,1:V)); 
[xc,Parac,Parac_cell,ffc,errc] = Cal_Loss(child_offspring);
child_data_Parac=Parac;
child_data_Parac_cell=Parac_cell;
error_norm=normalisation(errc);  %%                                
child_offspring=[xc ffc error_norm]; %% 

%% INtermediate population (Rt= Pt U Qt of 2N size)
data_Para_inter=[parent_data_Para; child_data_Parac]; %%
data_Para_cell_inter=[parent_data_Para_cell ; child_data_Parac_cell]; %%
population_inter=[population(:,1:V+M+1) ; child_offspring(:,1:V+M+1)]; %% 
% [population_inter_sorted,front]=NDS_CD_cons(population_inter);              % Non domination Sorting on offspring  
[population_inter_sorted,data_Para_sorted,data_Para_cell_sorted,front]=NDS_CD_cons(population_inter,data_Para_inter,data_Para_cell_inter);  
%% Replacement - N 

[new_pop,data_Para_new,data_Para_cell_new]=replacement(population_inter_sorted,data_Para_sorted,data_Para_cell_sorted, front);
population=real(new_pop);
parent_data_Para=data_Para_new;
parent_data_Para_cell=data_Para_cell_new;

end

[new_pop,sort]=sortrows(new_pop,V+1); 
Parar=parent_data_Para(sort,:); 
Parar_cell=parent_data_Para_cell(sort,:); 

Result_pop=new_pop(:,1:9);
[Result_population,Result_data_Parar_sorted,Result_data_Parar_cell_sorted,Result_front]=NDS_CD_cons(Result_pop, Parar, Parar_cell); 
Result_pop_feasible=Result_pop(find(population(:,V+M+1)==0),:);
Result_pop_infeasible=Result_pop(find(population(:,V+M+1)~=0),:);
flag=find(Result_pop(:,V+M+1)==0);


for i=1:size(Result_front(1).fr,2)

Front_pop(i,:)=Result_pop(flag(Result_front(1).fr(i)),:);
Front_Parac(i,:)=Parar(flag(Result_front(1).fr(i)),:); 
Front_Parac_cell(i,:)=Parar_cell(flag(Result_front(1).fr(i)),:); 
end
