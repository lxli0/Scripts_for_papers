%% The script estimates the temporal evolution of b-value for the catalog in Turkey
clc,clear
D=load('Ding_Turekey.txt');
Dt0=D(:,1);
jkf=find(Dt0>=2023.1013004791);
Dt=Dt0(jkf);
Magn0=D(jkf,6);
L0=length(Magn0);
bin=2000;
step=500;
num_bin=floor((L0-bin)/step)+1;
delta=0.05;
boostrap=100;
p=parpool(6);
parfor j=1:num_bin
    event_begin=(j-1)*step+1;
    event_end=(j-1)*step+bin;
    m_bin=Magn0(event_begin:event_end);
    t_bin=Dt(event_begin:event_end);
    
    [Mc_mle,sin_b_mle(j),number_mle]=MBS_MLE(m_bin);
    sigma_sin_mle(j)=sin_b_mle(j)/sqrt(number_mle);
    [Mc_pos,sin_b_pos(j),number_pos]=MBS_pos(m_bin);
    sigma_sin_pos(j)=sin_b_pos(j)/sqrt(number_pos);
    jkf1=find((m_bin-Mc_mle)>=-1e-6);
    m_bin_fil=m_bin(jkf1);
    sin_b_kms(j)=KMS(m_bin_fil); 
    sigma_sin_kms(j)=1.0297*sin_b_kms(j)/(number_mle)^0.4441;    
    
    b_mle0=zeros(1,boostrap);
    b_kms0=zeros(1,boostrap);
    b_pos0=zeros(1,boostrap);
    for k=1:boostrap
        random_Magn = m_bin(randi(length(m_bin),1,length(m_bin)));
        [Mc_mle,b_mle0(j,k),number_mle]=MBS_MLE(random_Magn);
        [Mc_pos,b_pos0(j,k),number_pos]=MBS_pos(random_Magn);
        jkf1=find((random_Magn-Mc_mle)>=-1e-6);
        random_Magn_fil=random_Magn(jkf1);
        b_kms0(j,k)=KMS(random_Magn_fil); 
    end
    mean_b_mle(j)=mean(b_mle0(j,:));
    std_b_mle(j)=std(b_mle0(j,:));
    mean_b_pos(j)=mean(b_pos0(j,:));
    std_b_pos(j)=std(b_pos0(j,:));
    mean_b_kms(j)=mean(b_kms0(j,:));
    std_b_kms(j)=std(b_kms0(j,:));

    t0(j)=t_bin(end);
    num_bin-j
end
delete(gcp('nocreate'))
t=(t0-2023.1013004791)*365.25;
%{
figure('units','normalized','position',[0.1,0.1,0.35,0.3])
color=1/255*[46 89 167;209 41 32;250 192 61];
h1=fill([t,fliplr(t)],[sin_b_mle+sigma_sin_mle,fliplr(sin_b_mle-sigma_sin_mle)],color(1,:));
set(h1,'edgealpha',0,'facealpha',0.3) 
hold on;
h11=plot(t,sin_b_mle,'color',color(1,:),'linewidth',2.5);
hold on;
h2=fill([t,fliplr(t)],[sin_b_pos+sigma_sin_pos,fliplr(sin_b_pos-sigma_sin_pos)],color(2,:));
set(h2,'edgealpha',0,'facealpha',0.3) 
hold on;
h22=plot(t,sin_b_pos,'color',color(2,:),'linewidth',2.5);
hold on;
h3=fill([t,fliplr(t)],[sin_b_kms+sigma_sin_kms,fliplr(sin_b_kms-sigma_sin_kms)],color(3,:));
set(h3,'edgealpha',0,'facealpha',0.3) 
hold on;
h33=plot(t,sin_b_kms,'color',color(3,:),'linewidth',2.5);
xlabel('Days Since 2023/2/7');
ylabel('b-value')
legend([h11,h22,h33],'MLE','b-positive','KMS','Location','South','NumColumns',3);
set(gca,'fontsize',16);
ylim([0.6 1.2])
grid on;
box on;
%}

figure('units','normalized','position',[0.1,0.1,0.5,0.3])
%color=1/255*[117 114 181; 197 86 89;203 180 123];
%color=1/255*[117 114 181;194 81 96;225 210 121];
color=1/255*[46 89 167;209 41 32;250 192 61];
h1=fill([t,fliplr(t)],[mean_b_mle+std_b_mle,fliplr(mean_b_mle-std_b_mle)],color(1,:));
set(h1,'edgealpha',0,'facealpha',0.3) 
hold on;
h11=plot(t,mean_b_mle,'color',color(1,:),'linewidth',2.5);
hold on;
h2=fill([t,fliplr(t)],[mean_b_pos+std_b_pos,fliplr(mean_b_pos-std_b_pos)],color(2,:));
set(h2,'edgealpha',0,'facealpha',0.3) 
hold on;
h22=plot(t,mean_b_pos,'color',color(2,:),'linewidth',2.5);
hold on;
h3=fill([t,fliplr(t)],[mean_b_kms+std_b_kms,fliplr(mean_b_kms-std_b_kms)],color(3,:));
set(h3,'edgealpha',0,'facealpha',0.3) 
hold on;
h33=plot(t,mean_b_kms,'color',color(3,:),'linewidth',2.5);
xlabel('Days Since 2023/2/7');
ylabel('b-value')
set(gca,'fontsize',16);
ylim([0.6 1.2])
grid on;
box on;
legend([h11,h22,h33],'MLE','b-positive','KMS','Location','South','NumColumns',3);

