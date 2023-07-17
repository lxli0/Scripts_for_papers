clc,clear
load('2to3_1delay.mat');% change for different GWS assumptions
A=load('Bao_TP_sensitivity.txt');% for the ETP
dh=100/1e3;
depth=dh/2:dh:100-dh/2;
depth0=A(:,2);
frequency=A(:,1);
sensitivity=A(:,3);
FRE=reshape(frequency,length(depth),[]);
SEN=reshape(sensitivity,length(depth),[]);
DEP=reshape(depth0,length(depth),[]);

[dt0,year2dt0,t0,F]=generate_waterload(350,1,0);

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
for j=1:length(Freq)
    [amplitude_loading(j),ind]=max(dc2c_loading(j,:));
    phase_loading(j)=t0(ind);
    depth_loading(j)=depth_loading0(j,ind);
    for k=1:3
        [amplitude_pp(j,k),ind]=max(dc2c_pp(j,:,k));
        phase_pp(j,k)=t0(ind);
        depth_pp(j,k)=depth_pp0(j,ind,k);
        [amplitude_eff(j,k),ind]=max(dc2c_eff(j,:,k));
        phase_eff(j,k)=t0(ind);
        depth_eff(j,k)=depth_eff0(j,ind,k);
    end
end

%% plot
%dc2c
figure('units','normalized','position',[0.1,0.1,0.6,0.2])
color=1/255*[
    192 137 44
    63 124 88
    30 87 181
    ];
subplot(1,2,1)
for j=1:3
    plot(Freq,phase_eff(:,j)/dt0/year2dt0*12,'-','color',color(j,:),'linewidth',2.5);
    hold on;
end
ylim([0 8])
set(gca,'ytick',0.5:1:7.5);
xlabel('Frequency (Hz)')
ylabel('Peak Time')
set(gca,'yticklabel',{'Jan',' ','Mar',' ','May',' ','Jul',' '}); 
box on;grid on;
set(gca,'fontsize',16);
xlim([min(Freq) max(Freq)])
%legend('D=0.1 m/s^2','D=1.0 m/s^2','D=10 m/s^2','Location','Southeast','NumColumns',2);

subplot(1,2,2)
for j=1:3
    plot(Freq,amplitude_eff(:,j),'-','color',color(j,:),'linewidth',2.5);
    hold on;
end
ylabel('Amplitude')
xlabel('Frequency (Hz)')
box on;grid on;
set(gca,'fontsize',16);
xlim([min(Freq) max(Freq)])
ylim([0 5e-3])

%origin depth
figure('units','normalized','position',[0.1,0.1,0.3,0.3])
for j=1:3
    plot(Freq,depth_eff(:,j),'-','color',color(j,:),'linewidth',2.5);
    hold on;
end
xlabel('Frequency (Hz)')
ylabel('Depth (km)')
box on;grid on;
set(gca,'fontsize',16);
xlim([min(Freq) max(Freq)])
ylim([0.6 1.6])
%legend('D=0.1 m/s^2','D=1.0 m/s^2','D=10 m/s^2');
