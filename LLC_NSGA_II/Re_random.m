function [x] = Re_random(x,xl,xu)
%RE_RANDOM 此处显示有关此函数的摘要
%   此处显示详细说明
global pop_size Vo

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
Flag_random=zeros(pop_size,1);
Vin=zeros(pop_size,1);
Vlm_pn=zeros(pop_size,1);
Random_times=zeros(pop_size,1);


for j=1:100
   
    n=x(:,1);
    Lr=x(:,2);
    fr=x(:,3);
    Lm=x(:,2).*x(:,4);%也就是x第四列是励磁电感与谐振电感比值k
    
    for i=1:pop_size
        Vin(i)=Vo/n(i);
        Cr(i)=1/4/pi/pi/fr(i)/fr(i)/Lr(i);
        Ln(i)=Lm(i)/Lr(i);
        Z0(i)=sqrt(Lr(i)/Cr(i));wr0(i)=1/sqrt(Lr(i)*Cr(i));
        Z1(i)=sqrt((Lr(i)+Lm(i))/Cr(i));wr1(i)=1/sqrt((Lr(i)+Lm(i))*Cr(i));
        wr(i)=wr0(i);
        wn(i)=wr1(i);
        In(i)=Vin(i)/Z0(i);
        Irm_estimate(i)=abs(Vo/Lm(i)/4/fr(i)/n(i));%以三角波来估计
        Ir0(i)=-Irm_estimate(i);
        Vcr0(i)=(Vin(i)-Vo/n(i)-n(i)*pi*1000*Z0(i)/2/Vo);%%vin-vo/n-npiPoZ0/2vo
        Ir0N(i)=Ir0(i)/In(i);
        Vcr0N(i)=Vcr0(i)/Vin(i);
        Ir2(i)=Irm_estimate(i);
        Vcr2(i)=-(Vin(i)-Vo/n(i)-n(i)*pi*1000*Z0(i)/2/Vo);
        Ir2N(i)=Ir2(i)/In(i);
        Vcr2N(i)=Vcr2(i)/Vin(i);
        t2(i)=0.5/fr(i);
        tr1(i)=1/fr(i)/2-t2(i);
        M(i)=1;
        r0(i)=sqrt(Ir0N(i)^2+(Vcr0N(i)-(1-M(i)))^2);
        the0(i)=atan(-Ir0N(i)/(Vcr0N(i)-(1-M(i))));
        r1(i)=sqrt((1+Ln(i))*Ir2N(i)^2+(Vcr2N(i)-1)^2);
        the1(i)=atan(-sqrt(1+Ln(i))*Ir2N(i)/(Vcr2N(i)-1));
        phi0(i)=t2(i)*wr0(i);
        phi1(i)=tr1(i)*wr1(i);  

        Vlm_pn(i)=abs(Lm(i)*(Vin(i)-Vin(i)*(r0(i)*cos(the0(i))+1-M(i)))/(Lr(i)+Lm(i)));%%为何电容电压前是-的符号，而且减去theta0
        if Vlm_pn(i)>0.8*Vo/n(i)
            Flag_random(i)=1;
        end
        if Flag_random(i)==1
            x(i,:) = xl+((xu-xl)*rand(1)); 
        end

    end
    Random_times(j)=sum(Flag_random);
    Flag_random=zeros(pop_size,1);
    if Random_times(j)==0
        break
    end
       
end
    
end

