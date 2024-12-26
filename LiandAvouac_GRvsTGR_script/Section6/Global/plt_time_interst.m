figure('units','normalized','position',[0.1,0.1,0.7,0.9])
time_node0=[1940 1970 2000 2013];
time_node=[min(T0),time_node0,max(T0)];

for i=1:length(time_node)-1
    jkf=find(T0>=time_node(i)&T0<time_node(i+1));
    M_sub=M0(jkf);
    subplot(3,length(time_node)-1,i)
    plt_duration_M(M_sub,Mc,delta)
    fortitle=[num2str(time_node(i)),' to ',num2str(time_node(i+1))];
    title(fortitle);
    subplot(3,length(time_node)-1,i+length(time_node)-1)
    plt_duration_DM(M_sub,delta)
end

function []=plt_duration_M(M_sub,Mc,delta)
    jkf=find(M_sub-Mc>=-1e-11);
    Dm_sub=M_sub(jkf);
    Dm_sub=Dm_sub-Mc;
    Mmin=0;
    [b_GR,loglikelihood_GR]=Estimation_GR_discrete(Dm_sub,Mmin,delta);
    [b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR_gridsearch_discrete(Dm_sub,Mmin,delta,b_GR);
    R=loglikelihood_TGR-loglikelihood_GR;
    P_LR=1-chi2cdf(2*R,1);
    plt_single_MFD(M_sub,Dm_sub,Mc,b_GR,b_TGR,mcorner_TGR,P_LR);
end



function []=plt_duration_DM(M_sub,delta)
    diffm=diff(M_sub);
    diffm0=diffm(diffm>=-1e-11);
    [Mc_prime,b_pos_1,number_abovemc]=MBS_MLE_discrete(diffm0,delta);
    diffm=diffm0(diffm0>=Mc_prime-1e-11);
    diffm=diffm-Mc_prime;
    [b_pos_2,loglikelihood_GR_dm]=Estimation_GR_discrete(diffm,0,delta);
    plt_single_MFD_DM(diffm0,diffm,Mc_prime,b_pos_2);
end
