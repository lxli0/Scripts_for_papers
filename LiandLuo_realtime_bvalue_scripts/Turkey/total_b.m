clc,clear
D=load('Ding_Turekey.txt');
Dt0=D(:,1);
jkf=find(Dt0>=2023.1013004791);
Dt=Dt0(jkf);
Dm=D(jkf,6);

[Mc_MBS,b_mle,number_mle]=MBS_MLE_discrete(Dm);
jkf=Dm>=(Mc_MBS-1e-6);
Dm_complete=Dm(jkf);
b_KMS=KMS(Dm_complete);
[Mc_pos,b_pos,number_pos]=MBS_pos_discrete(Dm);
