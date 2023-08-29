function [b_estimation_pos]=bpos_discrete(Dm0,Mmin)
    delta=0.05;
    b_estimation_pos=1/delta/log(10)*acoth((mean(Dm0)-Mmin+delta )/delta); 
end
