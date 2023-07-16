function [Tn,Magn]=gen_catalog(mp)
    dm=1.2;
    mmin=2;
    bd=1; % set b-value
    c=0.1; % day
    p=1.25;
    [NdT,md,td]=GR(bd,c,p,mmin,dm,mp,0);
    jkf=find(md>=mmin);
    se_num=0;
    while ~isempty(jkf)
        if se_num==0
            m_f=md(jkf);
            t_f=td(jkf);
        else
            m_f=md1(jkf);
            t_f=td1(jkf);
        end
        td1=[];md1=[];
        for j=1:length(jkf)
            md0=[];
            td0=[];
            [NdT,md0,td0]=GR(bd,c,p,mmin,dm,m_f(j),t_f(j));
            md1=[md1,md0];
            td1=[td1,td0];
        end
        md=[md,md1];
        td=[td,td1];
        jkf= find(md1>=mmin);
        se_num=se_num+1;
    end
    G=4.5;
    H=0.75;
    jkf0=find(td<=10); %filter the EQ by time
    T=td(jkf0);
    M=md(jkf0);
    T=[0,T];
    M=[mp,M];
    Matrix=[T;M]';
    sort_ma=sortrows(Matrix,1);
    T=sort_ma(:,1);
    M=sort_ma(:,2);
    % filter the earthquakes by STAI
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
    delta=0.05;
    Magn=roundn(Magn,log10(2*delta)); % magnitude discretization
    % generate first order aftershock for m0
    function [NdT,md,td]=GR(bd,c,p,mmin,dm,m0,t0)
      NdT=10^(bd*(m0-dm-mmin));
      md=zeros(1,floor(NdT));
        td=md;
        for i=1:floor(NdT)
            md(i)=mmin-log10(rand(1))/bd;
            td(i)=t0+(rand(1)^(1/(1-p))-1)*c;
        end
    end
end

