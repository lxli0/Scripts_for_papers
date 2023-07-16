clc,clear
figure('units','normalized','position',[0.1,0.1,0.45,0.35])
D=load('Bao_BB_velocity.txt');
D2=load('Bao_BB_sensitivity.txt');
load('lajolla.mat')
lajolla=flip(lajolla, 1);

depth=cumsum(D(:,1));
vs=D(:,3);
vp=D(:,2);

FRE=reshape(D2(:,1),1000,[]);
DEP=reshape(D2(:,2),1e3,[]);
phase_velocity=reshape(D2(:,3),1e3,[]);
phase_velocity(phase_velocity<=0)=1e-6;


subplot(1,2,1)
plot(vs,(depth-0.05),'linewidth',2.5,'color',lajolla(75,:));
hold on;
plot(vp,depth-0.05,'linewidth',2.5,'color',lajolla(200,:));
legend('V_S','V_P','Location','Southwest')
xlabel('V (km/s)')
ylabel('Depth (km)')
box on;
grid on;
set(gca,'fontsize',16);
set(gca,'YDir','reverse')  
set(gca,'xaxislocation','top')
xlim([1 8])
set(gca,'position',[0.12 0.2 0.25 0.6])

subplot(1,2,2)
sm=pcolor(FRE,DEP,phase_velocity);
set(sm,'FaceColor','interp','EdgeColor','none')
colormap(lajolla);
h=colorbar;
set(gca,'ColorScale','log')
set(get(h,'label'),'string','K_V [(100 m)^-^1]','FontSize',16);
xlabel('Frequency (Hz)');
ylabel('Depth (km)')
set(gca,'fontsize',16);set(gca,'YDir','reverse')  
ylim([0 40])
caxis([1e-6 1e-1]); 
set(gca,'position',[0.5 0.2 0.35 0.6])
