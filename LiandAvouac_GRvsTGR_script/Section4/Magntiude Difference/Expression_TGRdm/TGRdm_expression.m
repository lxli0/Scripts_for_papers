%clc,clear
%{
b0=[0.6 1 1.4];
L0=[1000 500 200];
zeta0=[0.8 0.7 0.6];
times=1000;
mcorner_TGR_est=zeros(length(b0),times);


for k=1:length(b0)
    b=b0(k);
    L=L0(k);
    zeta=zeta0(k);
    parpool(8);
    parfor i=1:times
        %% generate catalogs
        beta=2/3*b;
        mmin=0;
        Mmin=10^(1.5*(mmin+6.07));
        mcorner_TGR=zeta*log10(L)/b;
        M_corner=10^(1.5*(mcorner_TGR+6.07));

        U10=rand(1,L);
        U20=rand(1,L);
        U1=Mmin-M_corner*log(U10);
        U2=Mmin*U20.^(-1/beta);
        DM=min(U1,U2);
        Dm=2/3*log10(DM)-6.07;    
        
        dDm=diff(Dm);
        dDm_pos=dDm(dDm>=0);

        b_pos=1/log(10)/mean(dDm_pos);   

        %% Estimation
        [b_TGR(k,i),mcorner_TGR_est(k,i),loglikelihood_TGR]=Estimation_TGR(dDm,mmin);
    end
        delete(gcp('nocreate'))

end
save('solution_TGRdm.mat')
%}

%% clc,clear
load('solution_TGRdm.mat')
figure('units','normalized','position',[0.1,0.1,0.5,0.4])
colorm=jet(2*length(b0));
for k=1:length(b0)
    dDm_pos=[];
    L=L0(k);
    b=b0(k);
    zeta=zeta0(k);
    for jk=1:100
        beta=2/3*b;
        mmin=0;
        Mmin=10^(1.5*(mmin+6.07));
        mcorner_TGR=zeta*log10(L)/b;
        M_corner=10^(1.5*(mcorner_TGR+6.07));

        U10=rand(1,L);
        U20=rand(1,L);
        U1=Mmin-M_corner*log(U10);
        U2=Mmin*U20.^(-1/beta);
        DM=min(U1,U2);
        Dm=2/3*log10(DM)-6.07;    
        
        dDm=diff(Dm);
        dDm_pos=[dDm_pos,dDm(dDm>=0)];
    end
    plotgr(dDm_pos-0.05,colorm(2*k-1,:),2.5);
    hold on;
    m_rang=0:0.01:4;
    Mmin=10^(1.5*(mmin+6.07));
    M_range=10.^(1.5*(m_rang+6.07)); 
    mcorner_TGR_est_mean=mean(mcorner_TGR_est(k,:));
    M_corner=10^(1.5*(mcorner_TGR_est_mean+6.07));
    beta=2/3*mean(b_TGR(k,:));
    semilogy(m_rang,L/2*(Mmin./(M_range)).^beta.*exp((Mmin-M_range)/M_corner),'-','color',colorm(2*k,:),'linewidth',2.5);
     leg{2*k-1} = sprintf('synthetic catalogs ($b=%.1f,\\\\ n=%d\\\\,  m_{corner}=%.2f$)', b0(k), L0(k), mcorner_TGR);
     leg{2*k} = sprintf('fitting curve ($b^*=%.1f,\\\\ n^*=%d\\\\, m_{corner}^*=%.2f$)', mean(b_TGR(k,:)), L0(k)/2, mcorner_TGR_est_mean);
end
legend(leg, 'Interpreter', 'latex') % Specify 'Interpreter' as 'latex'
ylim([1,inf]);
grid on;box on; grid minor;
xlabel('Magnitude Difference');
ylabel('Culmulative Number');
set(gca,'fontsize',16);

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
    pl=semilogy(m,cn0/100,'-o','color',color,'linewidth',linewidth);
    hold on; 
end