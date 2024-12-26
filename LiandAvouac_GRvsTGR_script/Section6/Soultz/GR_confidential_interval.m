function [pl]=GR_confidential_interval(b,N,mc_set)
    color1=1/255*[132 94 194;178 91 0;0 139 200];
    n0=1:1:N;
    dm=0.01;
    m_range=0:dm:3*log10(N)/b;
    P_large=10.^(-b*m_range);
    alpha=90;
    for j=1:length(n0)
        n=n0(j);
        P=binopdf(n,N,P_large);
        [z0,h0]=max(P);
        meadian(j)=m_range(h0);
        TT=cumsum(P)/sum(P);
        [z1,h1]=min((abs(TT-(1-alpha/100)/2)));
        [z2,h2]=min((abs(TT-(1+alpha/100)/2)));
        small(j)=m_range(h1);
        large(j)=m_range(h2);
    end
    pl=semilogy(small+mc_set,n0,'--','color',color1(2,:),'linewidth',1.5);
    hold on;
    semilogy(large+mc_set,n0,'--','color',color1(2,:),'linewidth',1.5);
    hold on;
end