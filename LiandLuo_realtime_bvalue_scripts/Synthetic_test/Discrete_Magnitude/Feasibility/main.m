%% The script describes how nonuniform magnitude errors bias the b-value estimation for the three method (discrete magnitude)
clc,clear
Sim=500;%number of synthetic catalogs
for i=1:Sim
    [T,M]=gen_catalog(7);
    jkf=find(T>=1);
    T=T(jkf);
    Magn0=M(jkf);
    L=length(Magn0);
    Error=[0.25 0.05];
    error1=normrnd(0,Error(1),[1 L]);
    error2=normrnd(0,Error(2),[1 L]);
    Magn_Thr=3;
    jkf1=Magn0<Magn_Thr;
    jkf2=Magn0>=Magn_Thr;
    Magn=[Magn0(jkf1)+error1(jkf1),Magn0(jkf2)+error2(jkf2)];
    randIndex_A = randperm(L);
    Magn = Magn(randIndex_A);
   
    [Mc_mle(i),b_mle(i),number_mle(i)]=MBS_MLE(Magn);
    [Mc_pos(i),b_pos(i),number_pos(i)]=MBS_pos(Magn);  
    Mc_kms(i)=Mc_mle(i);
    number_kms(i)=length(Magn((Magn-Mc_kms(i))>=-1e-6));
    b_kms(i)=KMS(Magn,Mc_kms(i)); 
    aic_pos(i)=-4*number_pos(i)*log(2)+4*number_pos(i)*log(b_pos(i)+1) ...
        -2*number_pos(i)*(log(b_pos(i))+log(1))-2;
    aic_mle(i)=-4*number_mle(i)*log(2)+4*number_mle(i)*log(b_mle(i)+1) ...
        -2*number_mle(i)*(log(b_mle(i))+log(1))-2;
    aic_kms(i)=-4*number_kms(i)*log(2)+4*number_kms(i)*log(b_kms(i)+1) ...
        -2*number_kms(i)*(log(b_kms(i))+log(1))-2;
    p_pos(i)=exp(-aic_pos(i)/2-2);
    p_mle(i)=exp(-aic_mle(i)/2-2);
    p_kms(i)=exp(-aic_kms(i)/2-2);
    Sim-i
end

Matrix=[b_mle',b_pos',b_kms'];
Matrix(Matrix<=0.9)=0;
Matrix(Matrix>0.9&Matrix<1.1)=1;
Matrix(Matrix>=1.1)=0;
mle_ratio=length(find(Matrix(:,1)==1))/Sim;
pos_ratio=length(find(Matrix(:,2)==1))/Sim;
kms_ratio=length(find(Matrix(:,3)==1))/Sim;
sMa=sum(Matrix,2);
final=tabulate(sMa);


Matrix2=[p_mle',p_pos',p_kms'];
Matrix2(Matrix2<=5e-2)=0;
Matrix2(Matrix2>5e-2)=1;
mle_ratio_p=length(find(Matrix2(:,1)==1))/Sim;
pos_ratio_p=length(find(Matrix2(:,2)==1))/Sim;
kms_ratio_p=length(find(Matrix2(:,3)==1))/Sim;

sMa2=sum(Matrix2,2);
final2=tabulate(sMa2);

