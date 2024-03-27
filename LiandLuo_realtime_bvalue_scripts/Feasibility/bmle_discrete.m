function [b_estimation_MLE]=bmle_discrete(Dm0,Mmin)
    b_estimation_MLE=1/(log(10)*(mean(Dm0)-Mmin+0.05));
end
