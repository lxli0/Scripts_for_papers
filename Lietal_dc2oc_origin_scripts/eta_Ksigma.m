clc,clear
figure('units','normalized','position',[0.1,0.1,0.5,0.3])

A=load('Shen_BB_sensitivity.txt'); % change for different velocity profiles
depth0=A(:,2);
frequency=A(:,1);
sensitivity=A(:,3);

dh=100/1e3;
depth=dh/2:dh:100-dh/2;
rho=2700;
g=9.8;
lamda=0.37;
P0=(1-lamda)*g*rho;
P=depth*P0*1e3;
P=P/1e6;
logP=log10(P);
y1=10.^(-1.347*logP-6.882);
eta=y1';

FRE=reshape(frequency,length(depth),[]);
K_V=reshape(sensitivity,length(depth),[]);
DEP=reshape(depth0,length(depth),[]);

for i=1:length(frequency)/length(depth)
    K_sigma(:,i)=K_V(:,i).*eta;
end


subplot(1,2,1)
semilogx(eta,depth,'k','linewidth',2.5);
ylim([0 10])
set(gca,'YDir','reverse')  
set(gca,'xaxislocation','top')
grid on;
box on;
set(gca,'fontsize',16);
set(gca,'position',[0.12 0.2 0.25 0.6])
xlabel('\eta (Pa^-^1)')
xlim([1e-11 1e-5])
xticks([1e-11   1e-9 1e-7  1e-5])
ylabel('Depth (km)')


subplot(1,2,2)
K_sigma(K_sigma<=1e-10)=1e-10;
sm=pcolor(FRE,DEP,K_sigma/100);
set(sm,'FaceColor','interp','EdgeColor','none')
load('batlow.mat')
colormap(batlow);
h=colorbar;
set(gca,'ColorScale','log')
set(get(h,'label'),'string','K_\sigma (m^-^1Pa^-^1)','fontsize',16);
xlabel('Frequency (Hz)');
ylabel('Depth (km)')
caxis([1e-12 2e-11]); 
set(gca,'fontsize',16);
set(gca,'YDir','reverse') 
ylim([0 2])
set(gca,'position',[0.5 0.2 0.35 0.6])



