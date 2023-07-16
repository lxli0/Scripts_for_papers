function [b_kms] =KMS(Dm0,Mmin)
    Magn=Dm0((Dm0-Mmin)>=1e-6);
    size=length(Magn);
    kms0=-15.15*(log10(size))^(-2.14)+11.85;
    kms=zeros(1,10);
    for l=1:10
        dt_Poisson=exprnd(1,1,size-1);
        t_Poisson=cumsum(dt_Poisson);
        t_Poisson=[0,t_Poisson];
        randIndex_A = randperm(size);
        Magn_cal= Magn(randIndex_A);
        kms(l)=VGA(t_Poisson,Magn);
    end
    im_kms=mean(kms);
    b_kms=im_kms/kms0;
end