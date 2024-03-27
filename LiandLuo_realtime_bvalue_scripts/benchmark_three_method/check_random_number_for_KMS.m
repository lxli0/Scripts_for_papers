% The script tests the number of random catalog for KMS method
clc,clear
L=50; % number of events in each simulated catalog
Sim=1e4; % number of random catalogs
kms0=-15.15*(log10(L))^(-2.14)+11.85; % KMS/b ratio
b=1; % set b-value, we should expect the methods provide estimation close to this value
Magn_Start =0; % starting magnitude for simulations
b_es=zeros(1,Sim);
for i=1:Sim
    Pcm=rand(1,L);
    Magn=-1/b.*log10(Pcm)+ Magn_Start; % generate random magnitude
    kms=zeros(1,10);
    for j=1:10
        dt_Poisson=exprnd(1,1,L-1);
        t_Poisson=cumsum(dt_Poisson);
        t_Poisson=[0,t_Poisson]; % generate random occurrence time
        randIndex_A = randperm(L); 
        Magn_cal= Magn(randIndex_A); % randomly arrange the occurrence order
        kms(j)=VGA(t_Poisson,Magn); % calculate KMS for the generated catalog
        b_kms0(j)=kms(j)/kms0;
    end
    b_es0(i)=mean(b_kms0(1:9)); % estimated b-value by first 9 random catalogs
    b_es(i)=mean(b_kms0(1:10)); % estimated b-value by first 9 random catalogs
    per(i)=(b_es(i)-b_es0(i))/b_es0(i); %relative difference
    Sim-i
end
mean(abs(per)) %mean relative difference

% calculate KMS for a particular catalog
function [Slope] =VGA(Dt,Dm)
    L=length(Dm);
    VG=2*ones(1,L);
    VG(1)=1;
    VG(end)=1;
    for i=1:L-2
        for j=i+2:L
            tem_Dm=Dm(i+1:j-1);
            tem_Dt=Dt(i+1:j-1);
            cri_Dm=Dm(i)+(Dm(j)-Dm(i))*(tem_Dt-Dt(i))/(Dt(j)-Dt(i));
            cri=tem_Dm-cri_Dm;
            if max(cri)<0
                VG(i)=VG(i)+1;
                VG(j)=VG(j)+1;
            end
        end
    end
    p=polyfit(Dm,VG,1);
    Slope=p(1);
end
