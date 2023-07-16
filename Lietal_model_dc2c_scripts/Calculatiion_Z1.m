function [Z1]=Calculatiion_Z1(depth,eta)
    [dt0,year2dt0,t0,F]=generate_waterload(350,1,0);
    Z10=-2e-19;
    
    for j=1:length(depth)
        s11_0=-3/7;s22_0=-3/7;s33_0=-1;s12_0=0;s13_0=0;s23_0=0;% compression negative, niu=0.3
        s11=max(F)*s11_0;
        s22=max(F)*s22_0;
        s33=max(F)*s33_0;
        s12=max(F)*s12_0;
        s13=max(F)*s13_0;
        s23=max(F)*s23_0;
        delta_stress_3d_Boussinesq=[
           s11 s12 s13
           s12 s22 s23
           s13 s23 s33];
        [dsigma1,dsigma2,dsigma3]=principal_stress_change(delta_stress_3d_Boussinesq);
        [dVp1,dVp2,dVp3,dVs1,dVs2,dVs3]=stress2velocity(-dsigma1,-dsigma2,-dsigma3,Z10);
        dVs2Vs_loading(j)=(sqrt((2*(2e3+dVs3)^2+(2e3+dVs1)^2)/3)-2e3)/2e3;    
        Z1(j)=-s33*eta(j)/dVs2Vs_loading(j)*Z10;
    end
end