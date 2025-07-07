function [Flag_error,error_1,Flag_PN,Flag_OPO,Flag_PON,Flag_NOP,Flag_ZVS,Flag_M_less,Flag_M_more ]= NSGA_LLC_steady_cal_fs_all_range(n,Lr,fr,Lm,R)


global pop_size Vo Vin_max Vin_min  Coss tdead
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
x=zeros(7,pop_size);  


Flag_newton=zeros(pop_size,1);  %%% Newton iterative method error flag bit
Flag_M_less=zeros(pop_size,1);  %%% The flag bit that satisfies the gain M requirement of underresonance
Flag_M_more=zeros(pop_size,1);  %%% Marking bits that meet the gain M requirement of over-resonance
Flag_M=zeros(pop_size,1);       %%% The flag bit that meets the gain M requirement


Flag_PN=zeros(pop_size,1);
Flag_OPO=zeros(pop_size,1);
Flag_PON=zeros(pop_size,1);     
Flag_NOP=zeros(pop_size,1);
Flag_error=zeros(pop_size,1);   
Flag_ZVS =zeros(pop_size,1);       
error_1 =zeros(pop_size,1);   

ILs_rms_Matrix1_less    = zeros(pop_size,120);
Ipeak_t0t2_Matrix1_less = zeros(pop_size,120);
Ipeak_t2t3_Matrix1_less = zeros(pop_size,120);
Ipeak_Matrix1_less      = zeros(pop_size,120);
fs_Matrix1_less = zeros(pop_size,120);
r0_Matrix1_less = zeros(pop_size,120);
phi0_Matrix1_less   = zeros(pop_size,120);
the0_Matrix1_less   = zeros(pop_size,120);
r1_Matrix1_less = zeros(pop_size,120); 
phi1_Matrix1_less   = zeros(pop_size,120); 
the1_Matrix1_less   = zeros(pop_size,120);
M_Matrix1_less  = zeros(pop_size,120);
Vin_Matrix1_less    = zeros(pop_size,120);
ILs_rms_Matrix1_more    = zeros(pop_size,120);
Ipeak_Matrix1_more      = zeros(pop_size,120);
fs_Matrix1_more = zeros(pop_size,120);
r0_Matrix1_more = zeros(pop_size,120);
phi0_Matrix1_more   = zeros(pop_size,120);
the0_Matrix1_more   = zeros(pop_size,120);
r1_Matrix1_more     = zeros(pop_size,120);
phi1_Matrix1_more   = zeros(pop_size,120);
the1_Matrix1_more   = zeros(pop_size,120);
M_Matrix1_more  = zeros(pop_size,120);
Vin_Matrix1_more    = zeros(pop_size,120);

m=zeros(pop_size,1);
mn=zeros(pop_size,1);
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


for i=1:pop_size
    
   for j=1:1:10
    if ((n(i)*Vin_given(j))<=Vo)&&((n(i)*Vin_given(j+1))>Vo)
        mn(i)=j;
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

s1=0;
s2=0;


Vin(i)=Vo/n(i);
Cr(i)=1/4/pi/pi/fr(i)/fr(i)/Lr(i);
Ln(i)=Lm(i)/Lr(i);
Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i) );wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
wr(i)=wr0(i);
wn(i)=wr1(i);
In(i)=Vin(i)/Z0(i);

Flag_M_less(i)=0;  
Flag_M_more(i)=0;
Flag_M(i)=0;


R_OPO=n(i)*n(i)*Vo*pi*Z0(i)/2/(Vo*(Lm(i)+Lr(i))/Lm(i)-n(i)*Vin(i)) ;

if R>=R_OPO*0.95            
    Flag_OPO(i)=1;
    error_1(i)=10;
    Flag_error(i)=1;
    continue               
end

Vm=Lm(i)/(Lm(i)+Lr(i))*(2*R*Vin(i)-n(i)*Vo*pi*Z0(i))/2/R;
if Vm <= -Vo/n(i)
    Flag_PN(i)=1;
    error_1(i)=100;
    Flag_error(i)=1;
    continue                
end

%% newton iteration method

%%%%%% Before iteratively solving the steady-state value, the iterative initial value in the quasi-resonant state is estimated according to the transmitted parameters.  %%%%%%

Irm_estimate(i)=abs(Vo/Lm(i)/4/fr(i)/n(i));
Ir0(i)=-Irm_estimate(i);
Vcr0(i)=(Vin(i)-Vo/n(i)-n(i)*pi*1000*Z0(i)/2/Vo);
Ir0N(i)=Ir0(i)/In(i);
Vcr0N(i)=Vcr0(i)/Vin(i);
Ir2(i)=Irm_estimate(i);
Vcr2(i)=-(Vin(i)-Vo/n(i)-n(i)*pi*1000*Z0(i)/2/Vo);
Ir2N(i)=Ir2(i)/In(i);
Vcr2N(i)=Vcr2(i)/Vin(i);
t2(i)=0.5/fr(i);
tr1(i)=1/fr(i)/2-t2(i);


eps1=1e-7;
eps2=1e-7;
lemda=0.1;

%% fs<fr--newton iteration method

M(i)=1; 
r0(i)=sqrt(Ir0N(i)^2+(Vcr0N(i)-(1-M(i)))^2);
the0(i)=atan(-Ir0N(i)/(Vcr0N(i)-(1-M(i))));
r1(i)=sqrt((1+Ln(i))*Ir2N(i)^2+(Vcr2N(i)-1)^2);
the1(i)=atan(-sqrt(1+Ln(i))*Ir2N(i)/(Vcr2N(i)-1));
phi0(i)=t2(i)*wr0(i);
phi1(i)=tr1(i)*wr1(i);   

x(1:7,i)=[r0(i),phi0(i),the0(i),r1(i),phi1(i),the1(i),M(i)]; % 7xPopsize

step=-500;
j=1;

for j_fs=(fr(i)):step:0

    Ts(i,j)=1/j_fs;
    fs_Matrix1_less(i,j)=j_fs;
    iteration_times(i)=iteration_times(i)+1;
    if j_fs==(fr(i)-500)
        x(2,i)=3;
        x(5,i)=pi-3;
    end

        for k=1:1000
        
        b=F_less(x(:,i));

        norm_b=norm(b,inf);
        if norm_b<eps1 
            break;
        end
        A=Jac_less(x(:,i)); 
        dx=-pinv(A)*b; 
        x(:,i)=x(:,i)+lemda*dx; 

        norm_dx=norm(dx,inf);
        if norm_dx<eps2
            break;
        end
        end
     r0_Matrix1_less(i,j)=x(1,i);
    phi0_Matrix1_less(i,j)=x(2,i);
    the0_Matrix1_less(i,j)=x(3,i);
    r1_Matrix1_less(i,j)=x(4,i);
    phi1_Matrix1_less(i,j)=x(5,i);
    the1_Matrix1_less(i,j)=x(6,i);
    M_Matrix1_less(i,j)=x(7,i);    
    Vin_Matrix1_less(i,j)=Vo/n(i)/M_Matrix1_less(i,j);  

    if k==1000
            Flag_newton(i)=1; 
            error_1(i)=100;
            Flag_error(i)=1;
            break
    end

if abs(M_Matrix1_less(i,1)-1)>5e-3 || (r0_Matrix1_less(i,j)<0) || (r1_Matrix1_less(i,j)<0 )   
    error_1(i)=100;
    Flag_error(i)=1;
    break
end


    %%%%%% To determine whether the soft switching condition is satisfied %%%%%%
            if abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*sin(the0_Matrix1_less(i,j))/Z0(i))<(2*Vin_Matrix1_less(i,j)*Coss/tdead)
                Flag_ZVS(i)=1;
                error_1(i)=1;
                break
            end

    %%%%%%  Determine whether to enter the PON mode, that is, whether the current iteration result is reliable    %%%%%%
    
    %%%%% Calculate the minimum value of excitation inductance voltage %%%%%%%    
    if (the1_Matrix1_less(i,j)<pi && (phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))>pi)
        Vlm_min(i,j)=-Lm(i)*Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)/(Lm(i)+Lr(i));
    else
        Vlm_min(i,j)=min(Lm(i)*Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*cos(the1_Matrix1_less(i,j))/(Lm(i)+Lr(i)),Lm(i)*Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*cos(phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))/(Lm(i)+Lr(i)));  
    end 
    if Vlm_min(i,j)<=-Vo/n(i)

        Flag_PON(i)=1;
        error_1(i) =1;
        Flag_error(i)=1;
        break
    end


if (abs(Vin_given(m(i))-Vin_Matrix1_less(i,j))<2)&&(s1==0)
    r0_Matrix_less(i,m(i))=x(1,i);
    phi0_Matrix_less(i,m(i))=x(2,i);
    the0_Matrix_less(i,m(i))=x(3,i);
    r1_Matrix_less(i,m(i))=x(4,i);
    phi1_Matrix_less(i,m(i))=x(5,i);
    the1_Matrix_less(i,m(i))=x(6,i);
    M_Matrix_less(i,m(i))=x(7,i); 
    Vin_Matrix_less(i,m(i))=Vo/n(i)/M_Matrix_less(i,m(i));
    fs_Matrix_less(i,m(i))=fs_Matrix1_less(i,j);


  %%%%%Calculate the current value%%%%%%%        ------Including the average output current, the effective value of the resonant current and the peak value of the resonant current.

    ILs_rms_Matrix_less(i,m(i))=sqrt((r0_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*(2*phi0_Matrix1_less(i,j)+sin(2*the0_Matrix1_less(i,j))-sin(2*phi0_Matrix1_less(i,j)+2*the0_Matrix1_less(i,j)))/wr(i)/Ts(i,j)/Z0(i)/Z0(i)/2)+(r1_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*Vin_Matrix1_less(i,j)*(2*phi1_Matrix1_less(i,j)+sin(2*the1_Matrix1_less(i,j))-sin(2*phi1_Matrix1_less(i,j)+2*the1_Matrix1_less(i,j)))/wn(i)/Ts(i,j)/Z1(i)/Z1(i)/2));
    if (the0_Matrix1_less(i,j)<pi/2 && (phi0_Matrix1_less(i,j)+the0_Matrix1_less(i,j))>pi/2)||((the0_Matrix1_less(i,j)<-pi/2 && (phi0_Matrix1_less(i,j)+the0_Matrix1_less(i,j))>-pi/2))
        Ipeak_t0t2_Matrix_less(i,m(i))=abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)/Z0(i));
    else
        Ipeak_t0t2_Matrix_less(i,m(i))=max(abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*sin(the0_Matrix1_less(i,j))/Z0(i)),abs(Vin_Matrix1_less(i,j)*r0_Matrix1_less(i,j)*sin(phi0_Matrix1_less(i,j)+the0_Matrix1_less(i,j))/Z0(i)));
    end
    if (the1_Matrix1_less(i,j)<pi/2 && (phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))>pi/2)||((the1_Matrix1_less(i,j)<-pi/2 && (phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))>-pi/2))
        Ipeak_t2t3_Matrix_less(i,m(i))=abs(Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)/Z1(i));

    else
        Ipeak_t2t3_Matrix_less(i,m(i))=max(abs(Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*sin(the1_Matrix1_less(i,j))/Z1(i)),abs(Vin_Matrix1_less(i,j)*r1_Matrix1_less(i,j)*sin(phi1_Matrix1_less(i,j)+the1_Matrix1_less(i,j))/Z1(i)));

    end
    Ipeak_Matrix_less(i,m(i))=max(Ipeak_t0t2_Matrix_less(i,m(i)),Ipeak_t2t3_Matrix_less(i,m(i)));
    
    m(i)=m(i)-1;
    if(m(i)==0)
        m(i)=1;
        s1=1;
    end
end

    if M_Matrix1_less(i,j) > M_max(i)        
        Flag_M_less(i) = 1;
        break
    end

   j=j+1;
end


if Flag_error(i)==1 || Flag_M_less(i) == 0
    Flag_error(i)=1;
    continue
end

fs_min(i)=fs_Matrix1_less(i,j-1);



  if m(i)==11
        ILs_rms_max(i)=max(max(ILs_rms_Matrix1_less(i,:)),max(ILs_rms_Matrix1_more(i,:)));
        Ipeak_max(i)=max(max(Ipeak_Matrix1_less(i,:)),max(Ipeak_Matrix1_more(i,:)));
        continue
    end


%%   fs>fr--newton iteration method

M(i)=1;
r0(i)=sqrt(Ir0N(i)^2+(Vcr0N(i)-(1-M(i)))^2);
the0(i)=atan(-Ir0N(i)/(Vcr0N(i)-(1-M(i))));
r1(i)=sqrt(Ir2N(i)^2+(Vcr2N(i)+1+M(i))^2);
the1(i)=pi+atan(-Ir2N(i)/(Vcr2N(i)+1+M(i)));
phi0(i)=t2(i)*wr0(i);
phi1(i)=tr1(i)*wr0(i);

x(1:7,i)=[r0(i),phi0(i),the0(i),r1(i),phi1(i),the1(i),M(i)]; % 7xPopsize

j=1;

for j_fs=(fr(i)):500:2.5*fr(i)

    if Flag_error(i) == 1       
        break
    end 

    Ts(i,j)=1/j_fs;
    fs_Matrix1_more(i,j)=j_fs;
    iteration_times(i)=iteration_times(i)+1;
    
        for k=1:1000
        if k==1000
            Flag_newton(i)=1; 
            error_1(i) =100;
        break
        end

        b=F_more(x(:,i));
        norm_b=norm(b,inf);
        if norm_b<eps1 
            break;
        end
        A=Jac_more(x(:,i)); 
        dx=-pinv(A)*b; 
        x(:,i)=x(:,i)+lemda*dx; 
        norm_dx=norm(dx,inf);
        if norm_dx<eps2
            break;
        end
        end

        r0_Matrix1_more(i,j)=x(1,i);
        phi0_Matrix1_more(i,j)=x(2,i);
        the0_Matrix1_more(i,j)=x(3,i);
        r1_Matrix1_more(i,j)=x(4,i);
        phi1_Matrix1_more(i,j)=x(5,i);
        the1_Matrix1_more(i,j)=x(6,i);
        M_Matrix1_more(i,j)=x(7,i)    ;
        Vin_Matrix1_more(i,j)=Vo/n(i)/M_Matrix1_more(i,j); 



     
    
    %%%%%%  To determine whether the soft switching condition is satisfied.     %%%%%%
            if abs(Vin_Matrix1_more(i,j)*r0_Matrix1_more(i,j)*sin(the0_Matrix1_more(i,j)+phi0_Matrix1_more(i,j))/Z0(i))<(2*Vin_Matrix1_more(i,j)*Coss/tdead)
                Flag_ZVS(i)=1;
                error_1(i) =1;
                break
            end

        Vlm_N(i,j)=Lm(i)/(Lm(i)+Lr(i))*(Vin_Matrix1_more(i,j)-Vin_Matrix1_more(i,j)*(-r0_Matrix1_more(i,j)*cos(the0_Matrix1_more(i,j))+1-M_Matrix1_more(i,j)));
        if Vlm_N <=Vo/n(i)
            Flag_NOP(i) =1;
            Flag_error(i)=1;
            error_1(i) =1;
            break
        end

if(abs(Vin_given(mn(i)+1)-Vin_Matrix1_more(i,j))<2.5)&&(s2==0)
    r0_Matrix_more(i,mn(i)+1)=x(1,i);
    phi0_Matrix_more(i,mn(i)+1)=x(2,i);
    the0_Matrix_more(i,mn(i)+1)=x(3,i);
    r1_Matrix_more(i,mn(i)+1)=x(4,i);
    phi1_Matrix_more(i,mn(i)+1)=x(5,i);
    the1_Matrix_more(i,mn(i)+1)=x(6,i);
    M_Matrix_more(i,mn(i)+1)=x(7,i); 
    Vin_Matrix_more(i,mn(i)+1)=Vo/n(i)/M_Matrix_more(i,mn(i)+1);
    fs_Matrix_more(i,mn(i)+1)=fs_Matrix1_more(i,j);


%%%%%Calculate the current value%%%%%%%        ------Including the average output current, the effective value of the resonant current and the peak value of the resonant current.
    ILs_rms_Matrix_more(i,mn(i)+1)=sqrt((r0_Matrix1_more(i,j)*r0_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*(2*phi0_Matrix1_more(i,j)+sin(2*the0_Matrix1_more(i,j))-sin(2*phi0_Matrix1_more(i,j)+2*the0_Matrix1_more(i,j)))/wr(i)/Ts(i,j)/Z0(i)/Z0(i)/2)+(r1_Matrix1_more(i,j)*r1_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*Vin_Matrix1_more(i,j)*(2*phi1_Matrix1_more(i,j)+sin(2*the1_Matrix1_more(i,j))-sin(2*phi1_Matrix1_more(i,j)+2*the1_Matrix1_more(i,j)))/wr(i)/Ts(i,j)/Z0(i)/Z0(i)/2));
   
    if (phi0_Matrix1_more(i,j)+the0_Matrix1_more(i,j))>pi/2 && (the0_Matrix1_more(i,j)<pi/2 )
        Ipeak_Matrix_more(i,mn(i)+1)=abs(Vin_Matrix1_more(i,j)*r0_Matrix1_more(i,j)/Z0(i));
    else
        Ipeak_Matrix_more(i,mn(i)+1)=abs(Vin_Matrix1_more(i,j)*r0_Matrix1_more(i,j)*sin(phi0_Matrix1_more(i,j)+the0_Matrix1_more(i,j))/Z0(i));
    end

    mn(i)=mn(i)+1;
    if(mn(i)==11)
        mn(i)=10;
        s2=1;
    end
end
      
        if M_Matrix1_more(i,j) < M_min(i)         
        Flag_M_more(i) = 1;
        break
        end
        
        j=j+1;
end

if Flag_M_more(i) == 0
    error_1(i) =1;
    Flag_error(i)=1;
end

if Flag_error(i)==1 || Flag_M_more(i) == 0
    continue
end

 %%% Calculate the maximum and maximum peak value of the effective value of the resonant current in the whole working area.

    ILs_rms_max(i)=max(max(ILs_rms_Matrix_less(i,:)),max(ILs_rms_Matrix_more(i,:)));
    Ipeak_max(i)=max(max(Ipeak_Matrix_less(i,:)),max(Ipeak_Matrix_more(i,:)));

fs_max(i)=max( max(fs_Matrix_less(i,:)),max(fs_Matrix_more(i,:)));


end



















%% Time domain equations for PO mode 

function y=F_less(x)
y=zeros(size(x));
r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);

Ln(i)=Lm(i)/Lr(i);
Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i));wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
In(i)=Vin(i)/Z0(i);
IoN(i)=2*R*In(i)/(Ts(i,j)*wr0(i)*Vin(i)*n(i)*n(i));

y(1)=r0*sin(phi0+the0)-r1/(sqrt(1+Ln(i)))*sin(the1);
y(2)=-r0*cos(phi0+the0)-M+r1*cos(the1);
y(3)=r1/(sqrt(1+Ln(i)))*sin(phi1+the1)+r0*sin(the0);
y(4)=-r1*cos(phi1+the1)-r0*cos(the0)+2-M;
y(5)=r0*sin(phi0+the0)-r0*sin(the0)-M/Ln(i)*phi0;
y(6)=M-IoN(i)*(r0*cos(the0)-r0*cos(phi0+the0)-r0*sin(the0)*phi0-M/(2*Ln(i))*phi0^2);
y(7)=phi0/wr0(i)+phi1/wr1(i)-Ts(i,j)/2;

end

function Mx=Jac_less(x)

r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);


Ln(i)=Lm(i)/Lr(i);
Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i));wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
In(i)=Vin(i)/Z0(i);
IoN(i)=2*R*In(i)/(Ts(i,j)*wr0(i)*Vin(i)*n(i)*n(i));

Mx=[sin(phi0+the0) r0*cos(phi0+the0) r0*cos(phi0+the0) -1/sqrt(1+Ln(i))*sin(the1) 0 r1/sqrt(1+Ln(i))*cos(the1) 0;...
    -cos(phi0+the0) r0*sin(phi0+the0) r0*sin(phi0+the0) cos(the1) 0 -r1*sin(the1) -1;...
    sin(the0) 0 r0*cos(the0) 1/sqrt(1+Ln(i))*sin(phi1+the1) r1/sqrt(1+Ln(i))*cos(phi1+the1) r1/sqrt(1+Ln(i))*cos(phi1+the1) 0;...
    -cos(the0) 0 r0*sin(the0) -cos(phi1+the1) r1*sin(phi1+the1) r1*sin(phi1+the1) -1;...
    sin(phi0+the0)-sin(the0) r0*cos(phi0+the0)-M/Ln(i) r0*cos(phi0+the0)-r0*cos(the0) 0 0 0 -phi0/Ln(i);...
    -IoN(i)*(cos(the0)-cos(phi0+the0)-sin(the0)*phi0) -IoN(i)*(r0*sin(phi0+the0)-r0*sin(the0)-M/Ln(i)*phi0) -IoN(i)*(-r0*sin(the0)+r0*sin(phi0+the0)-r0*cos(the0)*phi0) 0 0 0 1-IoN(i)*(-1/(2*Ln(i))*phi0^2);...
    0 1/wr0(i) 0 0 1/wr1(i) 0 0]; 


end



%% Time domain equations for NP mode 

function y=F_more(x) 

y=zeros(size(x));
r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);

Ln(i)=Lm(i)/Lr(i);
Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i));wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
In(i)=Vin(i)/Z0(i);
IoN(i)=2*R*In(i)/(Ts(i,j)*wr0(i)*Vin(i)*n(i)*n(i));


y(1)=r0*sin(the0)+M*wr0(i)*Ts(i,j)/(4*Ln(i));
y(2)=r0*sin(phi0+the0)-r1*sin(the1);
y(3)=-r0*cos(phi0+the0)+r1*cos(the1)+2;
y(4)=r1*sin(phi1+the1)+r0*sin(the0);
y(5)=-r1*cos(phi1+the1)-r0*cos(the0)-2*M;
y(6)=M-IoN(i)*(r0*cos(the0)-r0*cos(phi0+the0)+r1*cos(the1)-r1*cos(phi1+the1));
y(7)=phi0/wr0(i)+phi1/wr0(i)-Ts(i,j)/2;


end


function Mx=Jac_more(x)

r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);

Ln(i)=Lm(i)/Lr(i);
Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i));wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
In(i)=Vin(i)/Z0(i);
IoN(i)=2*R*In(i)/(Ts(i,j)*wr0(i)*Vin(i)*n(i)*n(i));


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