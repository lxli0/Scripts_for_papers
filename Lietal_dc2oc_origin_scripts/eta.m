clc,clear
dh=100/1e3;
depth=dh/2:dh:100-dh/2;
rho=2700;
g=9.8;
lamda0=[0 0.37 0.9];
for i=1:3
    lamda=lamda0(i);
    P0=(1-lamda)*g*rho;
    P=depth*P0*1e3;
    P=P/1e6;
    logP=log10(P);
    y1=10.^(-1.347*logP-6.882);
    y2=10.^(-0.185*logP.^2-logP-6.261);
    y3=10.^(0.220*logP.^2-1.543*logP-7.938);
    semilogx(y1,depth,'k','linewidth',2.5);
    hold on;
end
set(gca,'YDir','reverse')  
set(gca,'xaxislocation','top')
grid on;
box on;
set(gca,'fontsize',16);
xlabel('\eta_S (Pa^-^1)')
ylabel('Depth (km)')
ylim([0 10]);
xlim([1e-11 1e-5])
xticks([1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5])
