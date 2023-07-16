function [Mc,b_es,number]=MBS_MLE(Dm)
    dm=0.1;
    D00=sort(Dm,'descend');
    m=0:dm:floor(D00(50)/dm)*dm;
    for k=1:length(m)
        Mmin=m(k);
        jkf=Dm>=(Mmin-1e-6);
        Dm0=Dm(jkf);
        b0(k)=bmle(Dm0,Mmin);
        sigma(k)=b0(k)/sqrt(length(Dm0));
    end
    Delta=0.5;
    number=Delta/dm+1;
    b0_mean=[];
    for k=1:length(m)-number+1
        b0_mean(k)=mean(b0(k:k+number-1));
    end
    b02=b0(1:length(m)-number+1);
    threshold=abs(b02-b0_mean)-0.5*sigma(1:length(m)-number+1);
    jkf=find(threshold<0);
    if ~isempty(jkf)
        Mc=m(jkf(1));
        jkff=find(Dm>=(Mc-1e-6));
        number=length(jkff);
        b_es=b02(jkf(1));
    else         
        [~,ind]=min(threshold);
        Mc=m(ind);
        jkff=find(Dm>=(Mc-1e-6));
        number=length(jkff);
        b_es=b02(ind);
    end
end