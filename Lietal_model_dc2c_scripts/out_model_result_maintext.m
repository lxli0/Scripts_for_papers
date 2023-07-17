clc,clear
load('2to3_1delay.mat');% change for different GWS assumptions
A=load('Bao_TP_sensitivity.txt');% for the ETP, and use 'Bao_BB_sensitivity.txt' for the BB
dh=100/1e3;
depth=dh/2:dh:100-dh/2;
depth0=A(:,2);
frequency=A(:,1);
sensitivity=A(:,3);
FRE=reshape(frequency,length(depth),[]);
SEN=reshape(sensitivity,length(depth),[]);
DEP=reshape(depth0,length(depth),[]);

[dt0,year2dt0,t0,F]=generate_waterload(350,1,0);%TWS

dc2c_loading0=zeros(length(frequency)/length(depth),length(t0),length(depth));
dc2c_pp0=zeros(length(frequency)/length(depth),length(t0),3,length(depth));
dc2c_eff0=zeros(length(frequency)/length(depth),length(t0),3,length(depth));

%% cumsum dc/c
for i=1:length(frequency)/length(depth)
    for mn=1:length(depth)
        for j=1:length(t0)
            dc2c_loading00(i,j,mn)=dVs2Vs_loading(j,mn)*SEN(mn,i);
            for k=1:3
                dc2c_pp00(i,j,k,mn)=dVs2Vs_pp(j,mn,k)*SEN(mn,i);
                dc2c_eff00(i,j,k,mn)=dVs2Vs_eff(j,mn,k)*SEN(mn,i);
            end
        end      
    end
end

for i=1:length(frequency)/length(depth)
    for j=1:length(t0)
        dc2c_loading0(i,j,:)=cumsum(dc2c_loading00(i,j,:));
        for k=1:3
           dc2c_pp0(i,j,k,:)=cumsum(dc2c_pp00(i,j,k,:));
           dc2c_eff0(i,j,k,:)=cumsum(dc2c_eff00(i,j,k,:));
        end
     end
end

%% dc/c       
for i=1:length(frequency)/length(depth)
    for j=1:length(t0)
        dc2c_loading(i,j)=dVs2Vs_loading(j,:)*SEN(:,i);
        for k=1:3
            dc2c_pp(i,j,k)=dVs2Vs_pp(j,:,k)*SEN(:,i);
            dc2c_eff(i,j,k)=dVs2Vs_eff(j,:,k)*SEN(:,i);
        end
    end
end

%% origin depth
thre=0.9;
for i=1:length(frequency)/length(depth)
    for mn=1:length(depth)
        for j=1:length(t0)
            ratio_loading(i,j,mn)=dc2c_loading0(i,j,mn)/dc2c_loading(i,j);
            for k=1:3
                ratio_pp(i,j,k,mn)=dc2c_pp0(i,j,k,mn)/dc2c_pp(i,j,k);
                ratio_eff(i,j,k,mn)=dc2c_eff0(i,j,k,mn)/dc2c_eff(i,j,k);
            end            
        end
    end
    
    for j=1:length(t0)
        if ~isnan(ratio_loading(i,j,1))
            jkf1=find(ratio_loading(i,j,:)>=thre);
            depth_loading0(i,j)=depth(jkf1(1));
            for k=1:3
                jkf2=find(ratio_pp(i,j,k,:)>=thre);
                depth_pp0(i,j,k)=depth(jkf2(1));
                jkf3=find(ratio_eff(i,j,k,:)>=thre);
                depth_eff0(i,j,k)=depth(jkf3(1));
            end
        end        
    end
end

Freq=FRE(1,:);
for j=1:length(Freq) % D=1m/s^2
    [amplitude_loading(j),ind]=max(dc2c_loading(j,:));
    phase_loading(j)=t0(ind);
    depth_loading(j)=depth_loading0(j,ind);
   
    [amplitude_pp(j),ind]=max(dc2c_pp(j,:,2));
    phase_pp(j)=t0(ind);
    depth_pp(j)=depth_pp0(j,ind,2);
    
    [amplitude_eff(j),ind]=max(dc2c_eff(j,:,2));
    phase_eff(j)=t0(ind);
    depth_eff(j)=depth_eff0(j,ind,2);
end

%% plot
% dc2c
figure('units','normalized','position',[0.1,0.1,0.6,0.6])
subplot(2,2,1)
plot(t0/dt0/year2dt0*12,F/9.8,'-k','linewidth',2.5);
hold on;
     [dt0_pp,year2dt0_pp,t0_pp,F_pp]=generate_waterload(2/3*350,1,1);
     plot(t0/dt0/year2dt0*12,F_pp/9.8,'--k','linewidth',2.5);
    hold on;

set(gca,'xtick',0.5:1:11.5);
set(gca,'xticklabel',{'Jan',' ','Mar',' ','May',' ','Jul',' ','Sep',' ','Nov',' '}); 
box on;
grid on;
set(gca,'fontsize',16);
ylim([-400,400])
legend('Terrestrial Water Storage','Groundwater Storage','NumColumns',1,'Location','South');
xlabel('Month')
xlim([0 12])
ylabel('EWH (mm)')

subplot(2,2,2)
Freq_plot=ones(length(t0),1)*Freq;
t0_plot=t0'*ones(1,length(Freq))/dt0/year2dt0*12;
sm=pcolor(t0_plot',Freq_plot',dc2c_loading);
set(sm,'FaceColor','interp','EdgeColor','none')
load('vik.mat');
colormap(vik);
%h=colorbar;
caxis([-4e-3 4e-3])
%set(get(h,'label'),'string','dc/c','FontSize',16);
set(gca,'xtick',0.5:1:11.5);
set(gca,'xticklabel',{'Jan',' ','Mar',' ','May',' ','Jul',' ','Sep',' ','Nov',' '}); 
xlabel('Month')
xlim([0 12])
ylabel('Frequency (Hz)')
set(gca,'fontsize',16);
grid on;
                
subplot(2,2,3)
sm=pcolor(t0_plot',Freq_plot',dc2c_pp(:,:,2));
set(sm,'FaceColor','interp','EdgeColor','none')
load('vik.mat');
colormap(vik);
%h=colorbar;
caxis([-4e-3 4e-3])
%set(get(h,'label'),'string','dc/c','FontSize',16);
set(gca,'xtick',0.5:1:11.5);
set(gca,'xticklabel',{'Jan',' ','Mar',' ','May',' ','Jul',' ','Sep',' ','Nov',' '}); 
xlabel('Month')
xlim([0 12])
ylabel('Frequency (Hz)')
set(gca,'fontsize',16);
grid on;

subplot(2,2,4)
sm=pcolor(t0_plot',Freq_plot',dc2c_eff(:,:,2));
set(sm,'FaceColor','interp','EdgeColor','none')
load('vik.mat');
colormap(vik);
%h=colorbar;
caxis([-4e-3 4e-3])
%set(get(h,'label'),'string','dc/c','FontSize',16);
set(gca,'xtick',0.5:1:11.5);
set(gca,'xticklabel',{'Jan',' ','Mar',' ','May',' ','Jul',' ','Sep',' ','Nov',' '}); 
xlabel('Month')
xlim([0 12])
ylabel('Frequency (Hz)')
set(gca,'fontsize',16);
grid on;

% origin depth
figure('units','normalized','position',[0.1,0.1,0.6,0.3])
color=1/255*[91 183 205; 197 86 89;203 180 123];
subplot(1,2,1)
plot(Freq,depth_loading,'color',color(1,:),'linewidth',2.5);
hold on
plot(Freq,depth_pp,'color',color(2,:),'linewidth',2.5);
hold on
plot(Freq,depth_eff,'color',color(3,:),'linewidth',2.5);
hold on
xlabel('Frequency (Hz)')
ylabel('Depth (km)')
box on;grid on;
set(gca,'fontsize',16);
xlim([min(Freq) max(Freq)])
legend('Elastic Load','Pore Pressure','Effective Stress','Location','Southeast')
ylim([0.6 1.6])
