function p=KStest(Dm,b)
    d=KSdistance(Dm,b);
    time=1000;
    d_theory=zeros(1,time);
    for i=1:time
        L=length(Dm);
        Pcm=rand(1,L);
        Mmin_true=0;
        Dm0=-1/b.*log10(Pcm)+Mmin_true;
        [b_GR,loglikelihood_GR]=Estimation_GR(Dm0,Mmin_true);
        d_theory(i)=KSdistance(Dm0,b_GR);
    end
    jkf=find(d_theory>=d);
    p=length(jkf)/time;
end

function d=KSdistance(Dm,b)
    L=length(Dm);
    dm=0.01;
    m=0:dm:max(Dm);
    n0=hist(Dm,m);
    cn0=zeros(1,length(m));
    cn0(1)=L;
    for i=2:length(m)
        cn0(i)=L-sum(n0(1:i-1)); 
    end
    Theoretical_cn0=L*10.^(-b*m);
    d=max(abs(cn0-Theoretical_cn0));
end