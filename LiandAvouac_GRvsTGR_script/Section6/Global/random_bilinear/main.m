clc,clear
close all
mag=load('Random_m.txt');
U=load('Random_U.txt');

figure('units','normalized','position',[0.1,0.1,0.7,0.25])
% plot theoretical distribution and random catalogs
subplot(1,3,1)
for i=1:100
Dm=gen_biGR(1e4);
plotgr(Dm,1/255*[255,222,173],1);
end
xlim([5 10])
ylim([1 1e4])
xlabel('Magnitude');
ylabel('Culmulative Number')
grid on;box on;grid minor;
set(gca,'fontsize',16)
hold on;
plot([5 7.5],[1e4 1e4*10^(-2.5)],'color',1/255*[255,165,0],'linewidth',2.5);
plot([7.5 10],[1e4*10^(-2.5) 1e4*10^(-2.5)*10^(-1.5*2.5)],'color',1/255*[255,165,0],'linewidth',2.5);

% plot two typical cases
subplot(1,3,2)
plt_detail_MFD(mag(:,1),5,0) % like TGR
subplot(1,3,3)
plt_detail_MFD(mag(:,2),5,0) % like GR

% check the random variable follows uniform distirbution
figure('units','normalized','position',[0.1,0.1,0.5,0.25])
subplot(1,2,1)
histogram(U(:,1))
xlabel('U')
ylabel('Count')
set(gca,'fontsize',16)
subplot(1,2,2)
histogram(U(:,2))
xlabel('U')
ylabel('Count')
set(gca,'fontsize',16)


% generate random samples from bi-linear GR distirbution
function [m]=gen_biGR(L)
    b1=1.0;
    b2=1.5;
    m=zeros(1,L);
    U0=rand(1,L);
    mc=5;
    deltam=0.05;
    mstart=mc-deltam; % for plotting purpose, otherwise, nonculmulative counts of m=6 will be smaller than expected
    mchange=7.6;
    m0=mchange-mc; %b1 for 
    threshold=10^(-m0*b1);
    for j=1:L
        if U0(j)>=threshold
            m(j)=-log10(U0(j))/b1;
        else
            m(j)=((b2-b1)*m0-log10(U0(j)))/b2;
        end
    end
    m=m+mstart;
end

% plot only culmulative MFD
function [pl]=plotgr(Dm,color,linewidth)
    Catalogsize=length(Dm);
    dm=0.1;
    m=0:dm:max(Dm);
    n0=hist(Dm,m); 
    cn0=[];
    cn0(1)=Catalogsize;
    for i=2:length(m)
        cn0(i)=Catalogsize-sum(n0(1:i-1)); 
    end
    pl=semilogy(m,cn0,'color',color,'linewidth',linewidth);
    hold on;
    
end

% plot detailed information about MFD
function []=plt_detail_MFD(M_sub,Mc,delta)
    jkf=find(M_sub-Mc>=-1e-11);
    Dm_sub=M_sub(jkf);
    Dm_sub=Dm_sub-Mc;
    Mmin=0;
    [b_GR,loglikelihood_GR]=Estimation_GR_discrete(Dm_sub,Mmin,delta);
    [b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR_gridsearch_discrete(Dm_sub,Mmin,delta,b_GR);
    R=loglikelihood_TGR-loglikelihood_GR;
    P_LR=1-chi2cdf(2*R,1);
    plt_single_MFD(M_sub,Dm_sub,Mc,b_GR,b_TGR,mcorner_TGR,P_LR);
end