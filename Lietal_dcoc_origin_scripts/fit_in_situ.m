clc,clear
h=[220 75 17 70 1e3 3.2e3 1e3];
eta=[1e-7 1e-6 1e-6 1e-7 2.4e-7 1e-9 2e-7];
lambda=[0 0.37 0.9];
for i=1:3
    P=2700*9.8*h;
    P=P*(1-lambda(i))/1e6; % in Mpa
    logP=log10(P);
    y1(i,:)=10.^(-1.347*logP-6.882);
    y2(i,:)=10.^(-0.185*logP.^2-logP-6.261);
    y3(i,:)=10.^(0.220*logP.^2-1.543*logP-7.938);
end
yy=[y1',y2',y3'];