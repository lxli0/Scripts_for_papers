function [dsigma1,dsigma2,dsigma3]=principal_stress_change(delta_stress)
    lij=[
        sqrt(3)/2 -1/2 0
        1/2  sqrt(3)/2  0
        0 0 1
        ];
    dsigma=lij'*delta_stress*lij;
    dsigma1=dsigma(1,1);
    dsigma2=dsigma(2,2);
    dsigma3=dsigma(3,3);
end
