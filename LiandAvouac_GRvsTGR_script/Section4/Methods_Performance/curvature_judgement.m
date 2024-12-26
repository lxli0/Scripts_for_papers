function p=curvature_judgement(Dm,b)
    c=curvature(Dm);
    time=1000;
    c_theory=zeros(1,time);
    for i=1:time
        L=length(Dm);
        Pcm=rand(1,L);
        Mmin_true=0;
        Dm0=-1/b.*log10(Pcm)+Mmin_true;
        c_theory(i)=curvature(Dm0);
    end
    jkf=find(c_theory>=c);
    p=length(jkf)/time;
end

function c=curvature(Dm)
    L=length(Dm);
    dm=0.01;
    m=0:dm:max(Dm);
    n0=hist(Dm,m);
    cn0=zeros(1,length(m));
    cn0(1)=L;
    for i=2:length(m)
        cn0(i)=L-sum(n0(1:i-1)); 
    end
    x0=m;
    y0=log10(cn0);
    p = polyfit(x0, y0, 2);
    x=linspace(min(x0),max(x0));
    y= polyval(p, x);
    x1 = diff(x);	
    x2 = diff(x1);	
    y1 = diff(y);
    y2 = diff(y1);
    x2(length(x1)) = x2(end);	
    y2(length(y1)) = y2(end);
    k = abs(x1.*y2-x2.*y1) ./ (x1.^2+y1.^2).^(3/2);
    k(length(x)) = k(end);
    c=mean(k);
end
