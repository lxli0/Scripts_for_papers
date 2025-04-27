clc,clear
figure('units','normalized','position',[0.1,0.1,0.2,0.6])
subplot(2,1,1)
D1=load('Bao_BB_input.txt');
depth=cumsum(D1(:,1));
vs=D1(:,3);
vp=D1(:,2);
plot(vs,depth-0.05,'linewidth',2.5,'color',1/255*[26 30 133]);
hold on;
plot(vp,depth-0.05,'linewidth',2.5,'color',1/255*[176 56 39]);
legend('V_S','V_P','Location','Southwest')
xlabel('V (km/s)')
ylabel('Depth (km)')
box on;grid on;
set(gca,'fontsize',16);
set(gca,'YDir','reverse')  
set(gca,'xaxislocation','top')
xlim([1 8])

subplot(2,1,2)
D1=load('Bao_TP_input.txt');
depth=cumsum(D1(:,1));
vs=D1(:,3);
vp=D1(:,2);
plot(vs,depth-0.05,'linewidth',2.5,'color',1/255*[26 30 133]);
hold on;
plot(vp,depth-0.05,'linewidth',2.5,'color',1/255*[176 56 39]);
legend('V_S','V_P','Location','Southwest')
xlabel('V (km/s)')
ylabel('Depth (km)')
box on;grid on;
set(gca,'fontsize',16);
set(gca,'YDir','reverse')  
set(gca,'xaxislocation','top')
xlim([1 8])