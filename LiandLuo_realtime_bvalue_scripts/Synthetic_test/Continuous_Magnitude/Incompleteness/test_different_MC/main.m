% The script evaluate the sensitivity of different methods on Mc/Mc' (continuous magnitudes)
clc,clear
figure('units','normalized','position',[0.1,0.1,0.4,0.3])
color=1/255*[46 89 167;209 41 32;250 192 61];
Sim = 50;% number of simulated catalogs
L = 3e3;% number of events in each simulated catalog
b=1;
Magn_Start = 0;
mu=0.8;
sigma0=0.2;
for i = 1 : Sim 
    b_kms=[];
    b_mle=[];
    b_pos=[];
    Pcm=rand(1,L);
    Magn0=-1/b.*log10(Pcm)+ Magn_Start;
    cri=normcdf(Magn0,mu,sigma0);
    ran_cri=rand(1,L);
    jkf=find(ran_cri<=cri);
    Magn=Magn0(jkf);
    dm=0.1;
    D00=sort(Magn,'descend');
    m1=0:dm:floor(D00(50)/dm)*dm;
    for k=1:length(m1)
        Mmin=m1(k);
        jkf=Magn>=(Mmin-1e-6);
        Dm0=Magn(jkf);
        b_mle(k)=1/(log(10)*(mean(Dm0)-Mmin));
        b_kms(k)=KMS(Dm0,Mmin);
    end
    %
    Dm=diff(Magn);
    D00=sort(Dm,'descend');
    m2=0:dm:floor(D00(50)/dm)*dm;
    for k=1:length(m2)
        Mmin=m2(k);
        jkf=Dm>=(Mmin-1e-6);
        Dm0=Dm(jkf);
        b_pos(k)=1/(log(10)*(mean(Dm0)-Mmin));
    end
    a=plot(m1,b_mle,'color',color(1,:),'linewidth',1.5);
   % a.Color(4)=0.5;
    hold on;
    a=plot(m2,b_pos,'color',color(2,:),'linewidth',1.5);
  %  a.Color(4)=0.5;
    hold on;
    a=plot(m1,b_kms,'color',color(3,:),'linewidth',1.5);
  %  a.Color(4)=0.5;
    hold on;
    legend('MLE','b-positive','KMS','Location','Southeast')
    Sim-i
end
xlabel('Mc or Mc''')
ylabel('b-value')
grid on;
box on;
xlim([0 2]);
ylim([0 2])
set(gca,'fontsize',16)
%{
subplot(1,2,1)
h1 = histogram( b_mle,20,'Normalization','probability'  ) ;
h1.EdgeColor=[1 1 1];
hold on;
h2 = histogram( b_pos,20,'Normalization','probability'  ) ;
h2.EdgeColor=[1 1 1];
h3 = histogram( b_kms,20,'Normalization','probability'  ) ;
h3.EdgeColor=[1 1 1];
xlabel( 'b-value' )
ylabel( 'Frequency' )
set(gca,'fontsize',16);
subplot(1,2,2)
h1 = histogram( Mc_mle,20,'Normalization','probability'  ) ;
h1.EdgeColor=[1 1 1];
hold on;
h2 = histogram( Mc_pos,20,'Normalization','probability'  ) ;
h2.EdgeColor=[1 1 1];
h3 = histogram( Mc_kms,20,'Normalization','probability'  ) ;
h3.EdgeColor=[1 1 1];
xlabel( 'Sample number' )
ylabel( 'Frequency' )
set(gca,'fontsize',16);
legend([h1,h2,h3],'MLE','b-positive','KMS');
grid on;
box on;
%nanmean(b_es)
%nanmean(Mc)
%}