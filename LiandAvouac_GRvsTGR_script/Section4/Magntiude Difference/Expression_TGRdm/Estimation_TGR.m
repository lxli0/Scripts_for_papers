function [b_TGR,Mcorner_TGR,loglikelihood_TGR]=Estimation_TGR(Dm,Mmin)
    DMoment=10.^(1.5*(Dm+6.07));
    Moment_min=10^(1.5*(Mmin+6.07));
    jkf=find((DMoment-Moment_min)>=-1e-16);
    n=length(jkf);
    
    A=1/n*sum(log(DMoment(jkf)/Moment_min));
    B=mean(DMoment(jkf))-Moment_min;

    syms x
    hx=1/n*sum(1./(1-x*(B-A.*DMoment(jkf))))-1;
    x0=NW(hx,1/B,100);

    Moment_corner_TGR=1/x0;
    Mcorner_TGR=log10(Moment_corner_TGR)*2/3-6.07;
    
    beta_TGR=(1-x0*B)/A;
    b_TGR=beta_TGR*3/2;
    loglikelihood_TGR=n*beta_TGR*log(Moment_min)+1/Moment_corner_TGR*(n*Moment_min-sum(DMoment(jkf))) ...
        -beta_TGR*sum(log(DMoment(jkf)))+sum(log(beta_TGR./DMoment(jkf)+1/Moment_corner_TGR));
end


function result=NW(h,x,n)
    f=matlabFunction(h); 
    f1=matlabFunction(diff(h));
    X(1)=x;
    i=2;
    while 1
        X(i)=X(i-1)-f(X(i-1))/f1(X(i-1));
        if abs(f(X(i))) <1e-5
             result=X(i);
             return;
        end
        if i>n
            result=X(i);
            return;
        end
        i=i+1;
    end
end