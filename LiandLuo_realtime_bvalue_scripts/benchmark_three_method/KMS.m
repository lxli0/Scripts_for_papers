function [b_kms]=KMS(Dm)
    Dm=reshape(Dm,1,[]);
    L=length(Dm);
    kms0=-15.15*(log10(L))^(-2.14)+11.85; % KMS/b ratio
    kms_every=zeros(1,10);
    for j=1:10
        dt_Poisson=exprnd(1,1,L-1);
        t_Poisson=cumsum(dt_Poisson);
        t_Poisson=[0,t_Poisson]; % generate random occurrence time
        randIndex_A = randperm(L); 
        Magn_cal= Dm(randIndex_A); % randomly arrange the occurrence order
        % calculate K for the generated catalog
        left=1;
        right=L;
        VG0=2*ones(1,L);
        VG0(1)=1;
        VG0(end)=1;
        VG=VGA(t_Poisson,Dm,left,right,VG0);
        % calculate KMS for the generated catalog
        pl=polyfit(Dm,VG,1);
        kms_every(j)=pl(1);
    end
    % calculate the corresponding b-value
    b_kms=mean(kms_every)/kms0;
end
