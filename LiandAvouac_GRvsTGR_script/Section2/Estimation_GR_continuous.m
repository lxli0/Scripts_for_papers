function [b_GR,loglikelihood_GR]=Estimation_GR_continuous(Dm,Mmin)
    DMoment=10.^(1.5*(Dm+6.07));
    Moment_min=10^(1.5*(Mmin+6.07));
    jkf=find((DMoment-Moment_min)>=-1e-16);
    n=length(jkf);
    beta_GR=n/sum(log( DMoment(jkf)/Moment_min ));
    b_GR=beta_GR*3/2;
    loglikelihood_GR=n*(beta_GR*log(Moment_min)+log(beta_GR))-(1+beta_GR)*sum(log(DMoment(jkf)));
end
