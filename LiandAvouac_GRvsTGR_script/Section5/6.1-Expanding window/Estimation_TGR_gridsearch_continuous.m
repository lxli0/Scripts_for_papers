function [b_TGR_output,Mcorner_TGR_output,loglikelihood_TGR_output]=Estimation_TGR_gridsearch_continuous(Dm,Mmin)
    DMoment=10.^(1.5*(Dm+6.07));
    Moment_min=10^(1.5*(Mmin+6.07));
    jkf=find((DMoment-Moment_min)>=-1e-16);
    n=length(jkf);
    Mcorner_TGR0=Mmin:0.05:ceil(max(Dm))+1;
    Moment_corner_TGR0=10.^(1.5*(Mcorner_TGR0+6.07));
    b_TGR0=0.6:0.05:1.5;    
    
    loglikelihood_TGR=1e9*ones(length(b_TGR0),length(Moment_corner_TGR0));
    for b_num=1:length(b_TGR0)
        b_TGR=b_TGR0(b_num);
        beta_TGR=2/3*b_TGR;
        for Mcorner_TGR_num=1:length(Moment_corner_TGR0)
            Moment_corner_TGR=Moment_corner_TGR0(Mcorner_TGR_num);
        loglikelihood_TGR(b_num,Mcorner_TGR_num)=n*beta_TGR*log(Moment_min)+1/Moment_corner_TGR*(n*Moment_min-sum(DMoment(jkf))) ...
            -beta_TGR*sum(log(DMoment(jkf)))+sum(log(beta_TGR./DMoment(jkf)+1/Moment_corner_TGR));  
        end
    end
    [minValue, linearIndex] = max(loglikelihood_TGR(:));
    [row, col] = ind2sub(size(loglikelihood_TGR), linearIndex);
    b_TGR_output=b_TGR0(row);
    Mcorner_TGR_output=Mcorner_TGR0(col);
    loglikelihood_TGR_output=minValue;
end

