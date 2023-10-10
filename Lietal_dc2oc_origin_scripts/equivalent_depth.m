figure('units','normalized','position',[0.1,0.1,0.9,0.6])
subplot(1,2,1)
x=5e-2:1e3;
x=1e6*x;
rho=2700;
g=9.8;
lamda=0.9;

xx=x/rho/g/(1-lamda);

loglog(xx,xx)
xlabel('Equivalent Depth (km)')
ylabel('\eta (Pa^-^1)')
grid on;
box on;
set(gca,'fontsize',16);

xlim([1e6*5e-2/rho/g/(1-lamda)/1e3,1e6*1e3/rho/g/(1-lamda)/1e3])
ylim([2e-12,inf])
set(gca,'ytick',[],'yticklabel',[]) 
