% This script provides b-value estimation for the three methods under different Mc/Mc' for the catalog in Italy 
%%
clc,clear
%% Load Data
Cat_Raw=load( 'Cat_16Nov_ZMAP.txt' ) ;
% select the events with the previous threshold criteria
Center   = [ 13.324 , 44.013 ] ;  % Long and Lat
DistMax  = 30;                   % Km
DepthMax = 30;                   % Km
Cat = Cat_Raw( Cat_Raw( : , 7 ) <= DepthMax & ...
      distance( Cat_Raw( : , 2 ) , Cat_Raw( : , 1 ) , Center(2) , Center(1) ) ...
      .* pi/180*6371 <= DistMax , : ) ;
% select events based on STAI
DeltaT = 1/6 ; %in days
Time = datenum( Cat( : , [ 3 : 5 , 8 : 10 ]) ) - ...
        datenum( Cat( 1 , [ 3 : 5 , 8 : 10 ]) ) ;
Cat_NoSTAI = Cat( Time >= DeltaT , : ) ;
Magn=Cat_NoSTAI(:,6);
%% plot magnitude-frequency distribution
Dm=Magn;
L=length(Dm);
dm=0.1;
m=0:dm:floor(max(Dm)/dm)*dm;
n0=hist(Dm,m);   
cn0(1)=L;
for i=2:floor(max(Dm)/dm)+1
    cn0(i)=L-sum(n0(1:i-1));
end
n=log10(n0);
cn=log10(cn0);
figure('units','normalized','position',[0.1,0.1,0.35,0.3])
color=1/255*[178 91 0;0 139 200];
plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 10, 'markerfacecolor', color(2,:), ...
            'MarkerEdgeColor', color(2,:));
hold on;
plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 10, 'markerfacecolor', color(1,:), ...
            'MarkerEdgeColor', color(1,:));
hold on;
legend([plot_n0,plot_cn0],'Non-cumulative','Cumulative','Location','Northeast');
grid on;
box on;
xlabel('Magnitude');
ylabel('Number');
set(gca,'fontsize',16);
xlim([0.9 inf]);
%% Calculation
delta=0.05;
D00=sort(Magn,'descend');
Mc_max=D00(50);
Mc_min=D00(end);
diff_M=diff(Magn);
D01=sort(diff_M,'descend');
Mc_prime_max=D01(50);
Mc_prime_min=0.1;
dm=0.1;
Mc=Mc_min:dm:Mc_max;
Mc_prime=Mc_prime_min:dm:Mc_prime_max;
boostrap=5000;
for i=1:length(Mc)
    jkf1=find((Magn-Mc(i))>=-1e-6);
    Magn=Magn(jkf1);
    %MLE
    sin_b_mle(i)=1/(log(10)*(mean(Magn)-Mc(i)+delta));
    sigma_sin_mle(i)=sin_b_mle(i)/sqrt(length(Magn));
    %KMS
    sin_b_kms(i)=KMS(Magn); 
    sigma_sin_kms(i)=1.0297*sin_b_kms(i)/(length(Magn))^0.4441;     
    
    b_mle0=zeros(1,boostrap);
    b_kms0=zeros(1,boostrap);
    for k=1:boostrap
        random_Magn = Magn(randi(length(Magn),1,length(Magn)));
        %MLE
        b_mle0(k)=1/(log(10)*(mean(random_Magn)-Mc(i)+delta));
        %KMS
        b_kms0(k)=KMS(random_Magn); 
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
    sin_b_pos(j)=1/delta/log(10)*acoth((mean(diff_M_pos)-Mc_prime(j)+delta )/delta); 
    sigma_sin_pos(j)=sin_b_pos(j)/sqrt(length(diff_M_pos));
    
    b_pos0=zeros(1,boostrap);
     for k=1:boostrap
        random_diff_M_pos = diff_M_pos(randi(length(diff_M_pos),1,length(diff_M_pos)));
        %b_positive
        b_pos0(k)=1/delta/log(10)*acoth((mean(random_diff_M_pos)-Mc_prime(j)+delta )/delta); 
     end
     mean_b_pos(j)=mean(b_pos0);
     std_b_pos(j)=std(b_pos0);
end
%% plot
%{
figure('units','normalized','position',[0.1,0.1,0.35,0.3])
color=1/255*[46 89 167;209 41 32;250 192 61];
h1=fill([Mc,fliplr(Mc)],[sin_b_mle+sigma_sin_mle,fliplr(sin_b_mle-sigma_sin_mle)],color(1,:));
set(h1,'edgealpha',0,'facealpha',0.3) 
hold on;
h11=plot(Mc,sin_b_mle,'color',color(1,:),'linewidth',2.5);
hold on;
h2=fill([Mc_prime,fliplr(Mc_prime)],[sin_b_pos+sigma_sin_pos,fliplr(sin_b_pos-sigma_sin_pos)],color(2,:));
set(h2,'edgealpha',0,'facealpha',0.3) 
hold on;
h22=plot(Mc_prime,sin_b_pos,'color',color(2,:),'linewidth',2.5);
hold on;
h3=fill([Mc,fliplr(Mc)],[sin_b_kms+sigma_sin_kms,fliplr(sin_b_kms-sigma_sin_kms)],color(3,:));
set(h3,'edgealpha',0,'facealpha',0.3) 
hold on;
h33=plot(Mc,sin_b_kms,'color',color(3,:),'linewidth',2.5);
xlabel('Mc or Mc''');
ylabel('b-value')
legend([h11,h22,h33],'MLE','b-positive','KMS','Location','Southeast');
set(gca,'fontsize',16);
ylim([0.2 1.2])
grid on;
box on;
%}

figure('units','normalized','position',[0.1,0.1,0.35,0.3])
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
set(gca,'fontsize',16);
grid on;
box on;
ylim([0.2 1.2])
legend([h11,h22,h33],'MLE','b-positive','KMS','Location','Southeast');
