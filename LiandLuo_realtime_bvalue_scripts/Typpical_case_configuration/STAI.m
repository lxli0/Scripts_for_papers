% Ths script provides an illustration of short-term aftershock incompleteness
clc,clear
figure('units','normalized','position',[0.1,0.1,0.25,0.25])

[T,M]=gen_catalog0(7);
scatter(1:length(M),M,20,'k','filled');
grid on;box on;
xlabel('Earthquake Number');ylabel('Magnitude');
set(gca,'fontsize',16);
xlim([0 1500])

function [Tn,Magn]=gen_catalog0(mp)
    m_min=2;
    [t_f,m_f]=BASS_primary(mp,0,m_min);
    jkf=m_f;
    md1=m_f;td1=t_f;
    while ~isempty(jkf)
        td=[];md=[];
        for j=1:length(jkf)
            [td0,md0]=BASS_primary(md1(j),td1(j),m_min);
            md=[md,md0];
            td=[td,td0];
        end
        md1=md;td1=td;
        m_f=[m_f,md];
        t_f=[t_f,td];
        jkf=md;
    end
   %% add STAI
    G=4.5;
    H=0.75;
    T=[0,t_f];
    M=[mp,m_f];
    Matrix=[T;M]';
    sort_ma=sortrows(Matrix,1);
    T=sort_ma(:,1);
    M=sort_ma(:,2);
    for h=2:length(M)-1
        miuM00=[];
        for j=1:h-1
            miuM00(j)=M(j)-G-H*log10(T(h)-T(j));
        end
        miuM0(h)=max(miuM00);
        miuM(h)=max(miuM0(h),2);
    end
    k=1;
    for h=2:length(M)-1
        cri=normcdf(M(h),miuM(h),0.2);
        ran_cri=rand(1);
        if ran_cri<=cri
            Magn(k)=M(h);
            Tn(k)=T(h);
            k=k+1;
        end
    end
    %% generate first order aftershock for m0
    function [t_primary,m_primary]=BASS_primary(m0,t0,m_min)
        dm=1.2;
        bd=1; % set b-value
        c=0.1; % day
        p=1.25;
        t_max=10; % day
        NdT=10^(bd*(m0-dm-m_min));
        U1=rand(1,floor(NdT));
        U2=rand(1,floor(NdT));
        m_primary=m_min-log10(U1)/bd;
        t_primary=t0+(U2.^(1/(1-p))-1)*c;
        jkf1=find(m_primary>=m_min);
        m_primary=m_primary(jkf1);
        t_primary=t_primary(jkf1);
        jkf2=find(t_primary<=t_max);
        m_primary=m_primary(jkf2);
        t_primary=t_primary(jkf2);
    end
end