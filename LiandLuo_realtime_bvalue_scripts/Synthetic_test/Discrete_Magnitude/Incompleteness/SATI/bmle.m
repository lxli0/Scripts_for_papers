function [b_estimation_MLE]=bmle(Dm0,Mmin)
    b_estimation_MLE=1/(log(10)*(mean(Dm0)-Mmin+0.05));
end
