clc,clear
color=1/255*[117 114 181;91 183 205; 197 86 89;203 180 123];
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
    P=depth*P0*1e3;%in Pa
    P=P/1e6;%in Mpa
    logP=log10(P);
    y1=10.^(-1.347*logP-6.882);
    y2=10.^(-0.185*logP.^2-logP-6.261);
    y3=10.^(0.220*logP.^2-1.543*logP-7.938);
    eta0(:,3*(i-1)+2)=y1';
    eta0(:,3*(i-1)+3)=y2';
    eta0(:,3*(i-1)+4)=y3';
end

%% KV
A1=load('Bao_TP_sensitivity.txt');
A2=load('Bao_BB_sensitivity.txt');
A3=load('Shen_TP_sensitivity.txt');
A4=load('Shen_BB_sensitivity.txt');

depth0=A1(:,2);
frequency=A1(:,1);
FRE=reshape(frequency,length(depth),[]);
DEP=reshape(depth0,length(depth),[]);

%% stress distribution
% uniform
dsigma1=1e3*ones(1,length(depth));

% EQ
Dsigma0=load('cross_out.txt');
Siz=size(Dsigma0);
kk=1001;% until 100 km
% center for fauly
jkf=find(Dsigma0(:,1)==1);
Dsigma=Dsigma0(jkf(1):jkf(1)+kk-1,1:6);
depth0=Dsigma(:,3);
dsigma20=-Dsigma(:,6);% compression is positive
dsigma2=interp1(depth0,1e5*dsigma20,depth,'spline');
%{
% A
Dsigma=Dsigma0(1:kk,1:6);
dsigma10=-Dsigma(:,6);% compression is positive
dsigma1=interp1(depth0,1e5*dsigma10,depth,'spline');

% B
Dsigma=Dsigma0(end-kk+1:end,1:6);
dsigma20=-Dsigma(:,6);% compression is positive
dsigma2=interp1(depth0,1e5*dsigma20,depth,'spline');
%}
% plot dsigma
figure('units','normalized','position',[0.05,0.1,0.7,0.35])
subplot(1,3,2)
plot(dsigma1,depth,'-','linewidth',2.5,'color',[0 0 0]);
grid on;box on;
set(gca,'YDir','reverse');
set(gca,'xaxislocation','top');
set(gca,'fontsize',16);
ylim([0 30])
ylabel('Depth (km)')
xlabel('\Delta \sigma_1 (Pa)')
set(gca,'position',[0.4 0.2 0.25 0.6])
subplot(1,3,3)
plot(dsigma2,depth,'-','linewidth',2.5,'color',[0 0 0]);
ylabel('Depth (km)');
xlabel('\Delta \sigma_2 (Pa)');
set(gca,'YDir','reverse');
set(gca,'xaxislocation','top');
set(gca,'fontsize',16);
grid on;box on;
ylim([0 30]);
set(gca,'position',[0.7 0.2 0.25 0.6]);

%% calculation for zh
ratio_thre=0.6;
for ijk=1:size(eta0,2)
    figure('units','normalized','position',[0.05,0.1,0.7,0.35])
    
    eta=eta0(:,ijk);
    
    for hh=1:4 % four velocity models
        if hh==1
            A=A1;
        elseif hh==2
            A=A2;
        elseif hh==3
            A=A3;
        else
            A=A4;
        end

        sensitivity=A(:,3);
        SEN=reshape(sensitivity,length(depth),[]);

        for i=1:length(frequency)/length(depth)
            Ksigma(:,i)=SEN(:,i).*eta;
            dc2c0_dsigma1(:,i)=Ksigma(:,i).*dsigma1';
            dc2c0_dsigma2(:,i)=Ksigma(:,i).*dsigma2';
        end
        R_dsigma1=cumsum(dc2c0_dsigma1,1)./sum(dc2c0_dsigma1,1);
        R_dsigma2=cumsum(dc2c0_dsigma2,1)./sum(dc2c0_dsigma2,1);
        for i=1:length(frequency)/length(depth)
            jkj1=find(R_dsigma1(:,i)>=ratio_thre);
            jkj2=find(R_dsigma2(:,i)>=ratio_thre);
            zh_dsigma1(i)=depth(jkj1(1));
            zh_dsigma2(i)=depth(jkj2(1));
        end

        subplot(1,3,1)
        semilogx(eta,depth,'k','linewidth',2.5);
        hold on;
        ylim([0 10])
        set(gca,'YDir','reverse')  
        set(gca,'xaxislocation','top')
        grid on;box on;
        set(gca,'fontsize',16);
        ylabel('Depth (km)')
        xlim([1e-11 1e-5])
        xticks([1e-11 1e-9 1e-7 1e-5])
        if ijk==1 || ijk==3
            xlabel('\eta_S (Pa^-^1)')
        else
            set(gca,'xticklabel',[]);
        end
        set(gca,'position',[0.1 0.2 0.25 0.6])    

        subplot(1,3,2)
        plot(FRE(1,:),zh_dsigma1,'-','linewidth',2.5,'color',color(hh,:));
        hold on;
        if hh==4
            xlim([0 1])
            xticks(0.1:0.2:0.9)
            ylim([0 1])
            ylabel('Depth (km)');
            if ijk==8 || ijk==10
                xlabel('Frequency (Hz)'); 
            else
                set(gca,'xticklabel',[]);
            end
            set(gca,'fontsize',16);
            grid on;box on;        
            set(gca,'position',[0.4 0.2 0.25 0.6])
        end
       
        subplot(1,3,3)
        plot(FRE(1,:),zh_dsigma2,'-','linewidth',2.5,'color',color(hh,:));
        hold on;  
        if hh==4
            xlim([0 1]) 
            xticks(0.1:0.2:0.9)
            ylim([0 1])
            ylabel('Depth (km)');
            if ijk==8 || ijk==10
                xlabel('Frequency (Hz)'); 
            else
                set(gca,'xticklabel',[]);
            end
            set(gca,'fontsize',16);
            grid on;box on;
            set(gca,'position',[0.7 0.2 0.25 0.6])
        end
    end
        
%    legend('ETP (Bao et al., 2015)','BB (Bao et al., 2015)', ...
%    'ETP (Shen et al., 2016)','BB (Shen et al., 2016)', ...
%   'Location','Southeast');
end