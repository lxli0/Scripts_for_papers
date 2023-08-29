function [b_estimation_MLE]=bmle_continuous(Dm0,Mmin)
    b_estimation_MLE=1/(log(10)*(mean(Dm0)-Mmin));
end
