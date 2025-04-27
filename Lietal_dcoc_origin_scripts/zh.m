clc,clear
close all
%% eta
dh=100/1e3;
depth=dh/2:dh:100-dh/2;
rho=2700;
g=9.8;

eta0(:,1)=1e-6*ones(length(depth),1);
lamda0=[0 0.37 0.9];
for i=1:3
    lamda=lamda0(i);
    P0=(1-lamda)*g*rho;
    P=depth*P0*1e3; % in Pa
    P=P/1e6; % in Mpa
    logP=log10(P);
    y1=10.^(-1.347*logP-6.882);
    y2=10.^(-0.185*logP.^2-logP-6.261);
    y3=10.^(0.220*logP.^2-1.543*logP-7.938);
    eta0(:,3*(i-1)+2)=y1';
    eta0(:,3*(i-1)+3)=y2';
    eta0(:,3*(i-1)+4)=y3';
    
end
eta0=[eta0(:,1:2),eta0(:,3),eta0(:,6),eta0(:,9),eta0(:,4),eta0(:,7),eta0(:,10)];
leg{1} = 'Uniform';
leg{2} = ['①, \lambda=', num2str(0)];
leg{3} = ['②, \lambda=', num2str(0)];
leg{4} = ['②, \lambda=', num2str(0.37)];
leg{5} = ['②, \lambda=', num2str(0.9)];
leg{6} = ['③, \lambda=', num2str(0)];
leg{7} = ['③, \lambda=', num2str(0.37)];
leg{8} = ['③, \lambda=', num2str(0.9)];

    
%% KV
A1=load('Bao_TP_sensitivity.txt');
A2=load('Bao_BB_sensitivity.txt');


depth0=A1(:,2);
frequency=A1(:,1);
FRE=reshape(frequency,length(depth),[]);
DEP=reshape(depth0,length(depth),[]);

%% stress distribution
% uniform
dsigma1=1e3*ones(1,length(depth));

% EQ related
Dsigma0=load('cross_out.txt');
Siz=size(Dsigma0);
kk=1001;% until 100 km

% center for fault
jkf=find(Dsigma0(:,1)==1);
Dsigma=Dsigma0(jkf(1):jkf(1)+kk-1,1:6);
depth0=Dsigma(:,3);
dsigma20=-Dsigma(:,6);% compression is positive
dsigma2=interp1(depth0,1e5*dsigma20,depth,'spline');

% A
Dsigma=Dsigma0(1:kk,1:6);
dsigma30=-Dsigma(:,6);% compression is positive
dsigma3=interp1(depth0,1e5*dsigma30,depth,'spline');

% B
Dsigma=Dsigma0(end-kk+1:end,1:6);
dsigma40=-Dsigma(:,6);% compression is positive
dsigma4=interp1(depth0,1e5*dsigma40,depth,'spline');


% plot dsigma
figure('units','normalized','position',[0.05,0.1,0.6,0.9])
subplot(3,3,2)
plot(dsigma1/1e3,depth,'-','linewidth',2.5,'color',[0 0 0]);
grid on;box on;
set(gca,'YDir','reverse');
set(gca,'xaxislocation','top');
set(gca,'fontsize',16);
xlim([0 2])
ylim([0 50])
ylabel('Depth (km)')
xlabel('$\Delta \sigma_1 $ (kPa)', 'Interpreter', 'latex');
subplot(3,3,3)
plot(dsigma2/1e6,depth,'-','linewidth',2.5,'color',[0 0 0]);
ylabel('Depth (km)');
xlabel('$\Delta \sigma_2 $ (MPa)', 'Interpreter', 'latex');
set(gca,'YDir','reverse');
set(gca,'xaxislocation','top');
set(gca,'fontsize',16);
grid on;box on;
ylim([0 50]);
xlim([-1 2])

%% calculation for z_h
ratio_thre=0.8;
color=1/255*[169,169,169
      255,215,0
    139, 0, 0;   % Dark Red
    205, 0, 0;   % Medium Red
    255, 99, 71; % Light Red
    0, 0, 139;   % Dark Blue
    	30,144,255;
    173, 216, 230 % Light Blue
    ];
for ijk=1:size(eta0,2)
    eta=eta0(:,ijk);
    
    for hh=1:2 % four velocity models
        if hh==1
            A=A1;
            D1=load('Bao_TP_input.txt'); % change for different profiles
        elseif hh==2
            A=A2;
        end

        sensitivity=A(:,3);
        K_V=reshape(sensitivity,length(depth),[]);

        for i=1:length(frequency)/length(depth)
            K_sigma(:,i)=K_V(:,i).*eta;
            dcoc0_dsigma1(:,i)=K_sigma(:,i).*dsigma1';
            dcoc0_dsigma2(:,i)=K_sigma(:,i).*dsigma2';
        end
        R_dsigma1=cumsum(dcoc0_dsigma1,1)./sum(dcoc0_dsigma1,1);
        R_dsigma2=cumsum(dcoc0_dsigma2,1)./sum(dcoc0_dsigma2,1);
        for i=1:length(frequency)/length(depth)
            jkj1=find(R_dsigma1(:,i)>=ratio_thre);
            jkj2=find(R_dsigma2(:,i)>=ratio_thre);
            z_h_dsigma1(i)=depth(jkj1(1));
            z_h_dsigma2(i)=depth(jkj2(1));
        end
        
        subplot(3,3,3+3*(hh-1)+1)
        FRE=reshape(A(:,1),1000,[]);
        DEP=reshape(A(:,2),1e3,[]);
        K_V=reshape(A(:,3),1e3,[]);
        K_V(K_V<=0)=1e-6;
        sm=pcolor(FRE,DEP,K_V/100);
        set(sm,'FaceColor','interp','EdgeColor','none')
        load('bamako.mat')
        colormap(bamako);
        h=colorbar;
        set(gca,'ColorScale','log')
        caxis([1e-8 1e-3])
        title(h, '$K_V$ (m$^{-1}$)', 'FontSize', 16, 'Interpreter', 'latex');
        xlabel('Frequency (Hz)');
        ylabel('Depth (km)')
        set(gca,'fontsize',16);set(gca,'YDir','reverse')  
        ylim([0 50])
        set(gca,'XTick', 0.2:0.2:0.8);
        set(gca,'xaxislocation','top')

        subplot(3,3,1)
        sw(ijk)=semilogx(eta,depth,'color',color(ijk,:),'linewidth',2.5);
        hold on;
        ylim([0 50])
        set(gca,'YDir','reverse')  
        set(gca,'xaxislocation','top')
        grid on;box on;grid minor;
        set(gca,'fontsize',16);
        ylabel('Depth (km)')
        xlim([1e-12 1e-5])
        xticks([1e-11 1e-9 1e-7 1e-5])
        xlabel('$\eta$ (Pa$^{-1}$)', 'Interpreter', 'latex');
        

        subplot(3,3,3+3*(hh-1)+2)
        semilogy(FRE(1,:), z_h_dsigma1, '-', 'linewidth', 2.5, 'color', color(ijk,:));
        hold on;
        xlim([0 1])
        ylim([1e-2 100])
        ylabel(['$z_h^{',num2str(100*ratio_thre),'\%}$ (km)'], 'Interpreter', 'latex');
        xlabel('Frequency (Hz)');
        set(gca, 'fontsize', 16);
        grid on; grid minor; box on;
        yticks([1e-2 1e-1 1e0 1e1 1e2])
        yticklabels({'10^{-2}', '10^{-1}', '10^{0}', '10^{1}', '10^{2}'})

        subplot(3,3,3+3*(hh-1)+3)
        semilogy(FRE(1,:), z_h_dsigma2, '-', 'linewidth', 2.5, 'color', color(ijk,:));
        hold on;  
        xlim([0 1])
        ylim([1e-2 100])
        ylabel(['$z_h^{',num2str(100*ratio_thre),'\%}$ (km)'], 'Interpreter', 'latex');
        xlabel('Frequency (Hz)');
        set(gca, 'fontsize', 16);
        grid on; grid minor; box on;
        yticks([1e-2 1e-1 1e0 1e1 1e2])
        yticklabels({'10^{-2}', '10^{-1}', '10^{0}', '10^{1}', '10^{2}'})
     
    end

end
subplot(3,3,1)
legend(sw,leg,'numcolumns',2,'location','south');