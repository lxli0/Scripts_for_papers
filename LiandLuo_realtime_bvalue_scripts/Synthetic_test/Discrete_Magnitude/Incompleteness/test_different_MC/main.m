%% The script evaluate the sensitivity of different methods on Mc/Mc' (discrete magnitudes)
clc,clear
figure('units','normalized','position',[0.1,0.1,0.4,0.3])
color=1/255*[46 89 167;209 41 32;250 192 61];
Sim = 50;%% number of simulated catalogs
L = 3e3;% number of events in each simulated catalog
b=1;
delta=0.05;
Magn_Start = -delta;
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
    Magn=roundn(Magn,log10(2*delta));
    dm=0.1;
    D00=sort(Magn,'descend');
    m1=0:dm:floor(D00(50)/dm)*dm;
    for k=1:length(m1)
        Mmin=m1(k);
        jkf=Magn>=(Mmin-1e-6);
        Dm0=Magn(jkf);
        b_mle(k)=1/(log(10)*(mean(Dm0)-Mmin+delta));
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
        b_pos(k)=1/delta/log(10)*acoth((mean(Dm0)-Mmin+delta )/delta);
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
ylim([0 2])
xlim([0 2])
grid on;
box on;
set(gca,'fontsize',16)