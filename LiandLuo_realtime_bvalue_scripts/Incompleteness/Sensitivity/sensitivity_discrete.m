%% The script evaluate the sensitivity of different methods on Mc/Mc' (discrete magnitudes)
clc,clear
figure('units','normalized','position',[0.1,0.1,0.9,0.3])
subplot(1,3,3);
color=1/255*[46 89 167;209 41 32;250 192 61];
Sim = 50;%% number of simulated catalogs
L = 3e3;% number of events in each simulated catalog
b=1;
delta=0.05;
Magn_Start = -delta;
mu=0.8;
sigma0=0.2;
Magn00=[];
Dm00=[];
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
    Magn00=[Magn00,Magn];
    dm=0.1;
    D00=sort(Magn,'descend');
    m1=0:dm:floor(D00(50)/dm)*dm;
    for k=1:length(m1)
        Mmin=m1(k);
        jkf=Magn>=(Mmin-1e-6);
        Dm0=Magn(jkf);
        b_mle(k)=1/(log(10)*(mean(Dm0)-Mmin+delta));
        b_kms(k)=KMS(Dm0);
    end
    %
    Dm=diff(Magn);
    Dm00=[Dm00,Dm];
    D00=sort(Dm,'descend');
    m2=0:dm:floor(D00(50)/dm)*dm;
    for k=1:length(m2)
        Mdmin=m2(k);
        jkf=Dm>=(Mdmin-1e-6);
        Dm0=Dm(jkf);
        b_pos(k)=1/delta/log(10)*acoth((mean(Dm0)-Mdmin+delta )/delta);
    end
    a1=plot(m1,b_mle,'color',color(1,:),'linewidth',1.5);
    a1.Color(4)=0.5;
    hold on;
    a2=plot(m2,b_pos,'color',color(2,:),'linewidth',1.5);
    a2.Color(4)=0.5;
    hold on;
    a3=plot(m1,b_kms,'color',color(3,:),'linewidth',1.5);
    a3.Color(4)=0.5;
    hold on;
    Sim-i
end
hold on;
plot( [ 0 2 ] , [ 1 1 ] , '--k' ,'linewidth',2.5 ) % add the real b-value line
legend([a1,a2,a3],'MLE','b-positive','KMS','Location','Northeast','NumColumns',3);
xlabel('Mc or Mc''')
ylabel('b-value')
grid on;
box on;
xlim([0 2]);
ylim([0 2])
set(gca,'fontsize',16)

%% Magnitude/Magnitude difference distribution 
subplot(1,3,1)
plot_dis(Magn00);
xlabel('Magnitude');
subplot(1,3,2)
dif_M=Dm00(Dm00>=-1e-6);
plot_dis(dif_M);
xlabel('Magnitude Difference');

%% function for plotting magnitude/magnitude difference  
function []=plot_dis(Dm)
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
    color=1/255*[178 91 0;0 139 200];
    plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
                'markersize', 7.5, 'markerfacecolor', color(2,:), ...
                'MarkerEdgeColor', color(2,:));
    hold on;
    plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
                'markersize', 7.5, 'markerfacecolor', color(1,:), ...
                'MarkerEdgeColor', color(1,:));
    hold on;
    legend([plot_n0,plot_cn0],'Non-cumulative','Cumulative','Location','Northeast');
    grid on;
    box on;
    ylabel('Number');
    set(gca,'fontsize',16);
    xlim([0 inf]);
end