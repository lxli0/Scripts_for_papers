clc,clear
close all

for i=1:10
    figure('units','normalized','position',[0.1,0.1,0.25,0.3])
    L=1e4;
    b1=1.0;
    b2=1.5;
    m=zeros(1,L);
    U0=rand(1,L);
    mc=5;
    deltam=0.05;
    mstart=mc-deltam; % for plotting purpose, otherwise, nonculmulative counts of m=6 will be smaller than expected
    mchange=7.6; %b1 for <mchange; b2 for >mchange
    m0=mchange-mc;
    threshold=10^(-m0*b1);
    for j=1:L
        if U0(j)>=threshold
            m(j)=-log10(U0(j))/b1;
        else
            m(j)=((b2-b1)*m0-log10(U0(j)))/b2;
        end
    end
    
    m=m+mstart;
    plt_detail_MFD(m,mc,0)
    hold on;
    xlim([5 10])
    m_store(i,:)=m;
    U0_stor(i,:)=U0;
end

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



