clc,clear
figure('units','normalized','position',[0.1,0.1,0.5,0.3])
D=load('Shen_BB_input.txt');
D2=load('Shen_BB_sensitivity.txt');


depth=cumsum(D(:,1));
vs=D(:,3);
vp=D(:,2);
subplot(1,2,1)
plot(vs,(depth-0.05),'linewidth',2.5,'color',1/255*[26 30 133]);
hold on;
plot(vp,depth-0.05,'linewidth',2.5,'color',1/255*[176 56 39]);
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
FRE=reshape(D2(:,1),1000,[]);
DEP=reshape(D2(:,2),1e3,[]);
phase_velocity=reshape(D2(:,3),1e3,[]);
phase_velocity(phase_velocity<=0)=1e-6;
sm=pcolor(FRE,DEP,phase_velocity/100);
set(sm,'FaceColor','interp','EdgeColor','none')
load('bamako.mat')
colormap(bamako);
h=colorbar;
set(gca,'ColorScale','log')
caxis([1e-8 1e-3])
set(get(h,'label'),'string','K_V (m^-^1)','FontSize',16);
xlabel('Frequency (Hz)');
ylabel('Depth (km)')
set(gca,'fontsize',16);set(gca,'YDir','reverse')  
ylim([0 40])
set(gca,'position',[0.5 0.2 0.35 0.6])

o1=DEP(:,1);
for i=1:20
    [val,pla]=max(phase_velocity(:,i));
    dep_max(i)=o1(pla);
end
