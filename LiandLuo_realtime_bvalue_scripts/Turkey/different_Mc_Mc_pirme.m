% This script provides b-value estimation for the three methods under different Mc/Mc' for the catalog in Italy 
%% Load Data
clc,clear
D=load('Ding_Turkey.txt');
Dt0=D(:,1);
jkf=find(Dt0>=2023.1013004791);
Dt=Dt0(jkf);
Magn=D(jkf,6);

%% Calculation
delta=0.05;
D00=sort(Magn,'descend');
Mc_max=D00(50);
Mc_min=0;
diff_M=diff(Magn);
D01=sort(diff_M,'descend');
Mc_prime_max=D01(50);
Mc_prime_min=0.1;
dm=0.1;
Mc=Mc_min:dm:Mc_max;
Mc_prime=Mc_prime_min:dm:Mc_prime_max;
for i=1:length(Mc)
    jkf1=find((Magn-Mc(i))>=-1e-6);
    Magn=Magn(jkf1);
    boostrap=floor(1e5/length(Magn));
    b_mle0=zeros(1,boostrap);
    b_kms0=zeros(1,boostrap);
    for k=1:boostrap
        random_Magn = Magn(randi(length(Magn),1,length(Magn)));
        %MLE
        b_mle0(k)=1/(log(10)*(mean(random_Magn)-Mc(i)+delta));
        %KMS
        b_kms0(k)=KMS(random_Magn); 
        boostrap-k
    end
    mean_b_mle(i)=mean(b_mle0);
    std_b_mle(i)=std(b_mle0);
    mean_b_kms(i)=mean(b_kms0);
    std_b_kms(i)=std(b_kms0);
end
for j=1:length(Mc_prime)
    jkf2=find((diff_M-Mc_prime(j))>=-1e-6);
    diff_M_pos=diff_M(jkf2);
    %b_positive
    b_pos0=zeros(1,boostrap);
     for k=1:boostrap
        random_diff_M_pos = diff_M_pos(randi(length(diff_M_pos),1,length(diff_M_pos)));
        b_pos0(k)=1/delta/log(10)*acoth((mean(random_diff_M_pos)-Mc_prime(j)+delta )/delta); 
     end
     mean_b_pos(j)=mean(b_pos0);
     std_b_pos(j)=std(b_pos0);
end
%% plot
figure('units','normalized','position',[0.1,0.1,0.5,0.3])
color=1/255*[46 89 167;209 41 32;250 192 61];
h1=fill([Mc,fliplr(Mc)],[mean_b_mle+std_b_mle,fliplr(mean_b_mle-std_b_mle)],color(1,:));
set(h1,'edgealpha',0,'facealpha',0.3) 
hold on;
h11=plot(Mc,mean_b_mle,'color',color(1,:),'linewidth',2.5);
hold on;
h2=fill([Mc_prime,fliplr(Mc_prime)],[mean_b_pos+std_b_pos,fliplr(mean_b_pos-std_b_pos)],color(2,:));
set(h2,'edgealpha',0,'facealpha',0.3) 
hold on;
h22=plot(Mc_prime,mean_b_pos,'color',color(2,:),'linewidth',2.5);
hold on;
h3=fill([Mc,fliplr(Mc)],[mean_b_kms+std_b_kms,fliplr(mean_b_kms-std_b_kms)],color(3,:));
set(h3,'edgealpha',0,'facealpha',0.3) 
hold on;
h33=plot(Mc,mean_b_kms,'color',color(3,:),'linewidth',2.5);
xlabel('Mc or Mc''');
ylabel('b-value')
set(gca,'fontsize',20);
grid on;
box on;
ylim([0.2 1.8])
legend([h11,h22,h33],'MLE','b-positive','KMS','Location','Southeast');
