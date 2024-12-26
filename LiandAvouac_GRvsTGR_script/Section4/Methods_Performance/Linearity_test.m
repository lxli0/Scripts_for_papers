function [judgement]=Linearity_test(Dm)
    dm=0.1;
    D00=sort(Dm,'descend');
    m=0:dm:floor(D00(50)/dm)*dm;
    for k=1:length(m)
        Mmin=m(k);
        [b_GR(k),loglikelihood_GR]=Estimation_GR(Dm,Mmin);
    end
    n=length(find((Dm-m(end))>=-1e-6));
    sigma_50=b_GR(end)/sqrt(n);
    NLINDEX=std(b_GR)/sigma_50;
    if NLINDEX<=1
        judgement=1; % GR
    else
        judgement=2; % TGR
    end
end