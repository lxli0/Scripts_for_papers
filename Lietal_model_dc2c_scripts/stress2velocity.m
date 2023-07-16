function [dVp1,dVp2,dVp3,dVs1,dVs2,dVs3] = stress2velocity(dsigma1,dsigma2,dsigma3,Z1)
    %% initial parameters
    rho=2400;
    Vp0=3000;
    Vs0=2000;
    %% C0 and S0
    C0=zeros(6,6);
    C0(1,1)=rho*Vp0^2;
    C0(2,2)=C0(1,1);
    C0(3,3)=C0(1,1);
    C0(1,2)=rho*Vp0^2-2*rho*Vs0^2;
    C0(1,3)=C0(1,2);
    C0(2,3)=C0(1,2);
    C0(4,4)=rho*Vs0^2;
    C0(5,5)=C0(4,4);
    C0(6,6)=C0(4,4);

    S0=inv(C0);

    %% final C and S

    daerfa11=2*pi/15*(3*dsigma1+dsigma2+dsigma3)*Z1;
    daerfa22=2*pi/15*(dsigma1+3*dsigma2+dsigma3)*Z1;
    daerfa33=2*pi/15*(dsigma1+dsigma2+3*dsigma3)*Z1;

    dS=zeros(6,6);
    dS(1,1)=daerfa11;
    dS(2,2)=daerfa22;
    dS(3,3)=daerfa33;
    dS(4,4)=daerfa22+daerfa33;
    dS(5,5)=daerfa11+daerfa33;
    dS(6,6)=daerfa11+daerfa22;

    S=dS+S0;
    C=inv(S);
    
    %% final velocity
    V11=sqrt(C(1,1)/rho);
    V12=sqrt(C(5,5)/rho);
    V13=sqrt(C(6,6)/rho);
    V21=sqrt(C(4,4)/rho);
    V22=sqrt(C(2,2)/rho);
    V23=sqrt(C(6,6)/rho);
    V31=sqrt(C(4,4)/rho);
    V32=sqrt(C(5,5)/rho);
    V33=sqrt(C(3,3)/rho);
    
    %% velocity changes
    dVp1=V11-Vp0;
    dVp2=V22-Vp0;
    dVp3=V33-Vp0;
    dVs1=V21-Vs0;
    dVs2=V12-Vs0;
    dVs3=V13-Vs0;
end