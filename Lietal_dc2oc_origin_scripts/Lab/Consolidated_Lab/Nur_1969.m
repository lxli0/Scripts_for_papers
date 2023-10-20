clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=[0
50
100
200
400
700
1000
2000
3000];
Velocity=[3.3000 	5.3000 	2.3200 	2.4200 	3.8000 	5.4800 	2.8000 	3.0000 	4.5000 	5.7000 	2.9000 	2.9000 	5.0000 	6.4000 	3.7000 	3.7000 	2.6000 	4.5400 	1.4400 	1.4600 	5.5900 	5.6200 	2.9900 	2.9800 
4.2000 	5.8100 	2.5600 	2.7600 	4.5000 	5.6100 	3.0300 	3.0700 	5.6500 	6.2200 	3.2300 	3.2300 	5.9000 	6.6000 	3.7500 	3.7300 	2.7900 	4.6400 	1.6400 	1.6800 	5.6100 	5.6300 	3.0000 	2.9900 
5.0500 	6.0200 	2.7900 	3.0000 	4.9800 	5.7000 	3.0700 	3.1000 	5.9100 	6.2200 	3.3300 	3.3300 	6.4000 	6.7300 	3.8400 	3.8300 	3.0400 	4.6800 	1.7500 	1.7600 	5.6200 	5.6400 	3.0100 	3.0100 
5.6900 	6.1400 	3.1200 	3.1600 	5.3100 	5.7700 	3.1700 	3.1500 	6.1700 	6.3400 	3.4200 	3.4200 	6.6900 	6.8100 	3.9700 	3.9300 	3.3500 	4.7600 	1.9200 	1.9400 	5.6300 	5.6700 	3.0300 	3.0400 
6.0000 	6.2700 	3.4400 	3.3200 	5.5500 	5.9000 	3.2800 	3.2300 	6.3400 	6.3800 	3.5100 	3.5100 	6.8700 	6.8900 	4.0600 	4.0700 	3.7700 	4.8200 	2.1400 	2.1400 	5.6400 	5.7000 	3.0500 	3.0600 
6.3000 	6.4400 	3.6000 	3.4800 	5.7300 	6.0100 	3.3500 	3.3100 	6.4300 	6.4500 	3.5600 	3.5600 	6.9200 	6.9600 	4.1000 	4.1500 	4.1000 	4.8600 	2.3200 	2.2800 	5.6600 	5.7300 	3.0800 	3.0800 
6.4600 	6.4800 	3.6600 	3.6000 	5.8500 	6.0600 	3.4000 	3.3500 	6.4500 	6.5100 	3.5800 	3.5800 	6.9400 	6.9900 	4.1200 	4.1900 	4.3800 	4.9000 	2.4700 	2.4000 	5.6800 	5.7500 	3.0900 	3.0900 
6.5500 	6.5400 	3.7300 	3.6900 	6.0600 	6.1300 	3.4800 	3.4400 	nan	nan	nan	nan	7.0000 	7.0500 	4.1600 	4.2600 	4.8300 	4.9600 	2.6800 	2.6100 	5.7200 	5.7800 	3.1000 	3.1000 
6.5800 	6.5800 	3.7600 	3.7300 	6.1300 	6.1700 	3.5200 	3.4700 	nan	nan	nan	nan	7.0300 	7.0900 	4.1900 	4.2900 	nan	nan	nan	nan	5.7500 	5.8000 	3.1200 	3.1100 
																							
];
Sample_number=size(Velocity,2);

%% create colormap
range=0:size(lajolla,1)-1;
range=range*(Sample_number-1)/(size(lajolla,1)-1)+1;
color=nan(Sample_number,3);
for i=1:3
    color(:,i)=interp1(range,lajolla(:,i),1:Sample_number);
end
%% Process
for i=1:Sample_number
    P0=1e5*Pressure;
    V0=Velocity(:,i);
    [mean_P00,eta00]=process(P0,V0);
    mean_P0(i,1:length(mean_P00))=mean_P00;
    eta0(i,1:length(mean_P00))=eta00;
    figure(2)
    semilogy(mean_P00,eta00,'o','markersize',8, 'markerfacecolor',color(i,:),'MarkerEdgeColor', color(i,:));
    hold on;
end
figure(1)
grid on;
box on;
set(gca,'fontsize',16);
xlabel('Pressure (Pa)')
ylabel('Velocity (km/s)')
figure(2)
eta=reshape(eta0,1,[])';
jk=find(eta>0);
eta=eta(jk);
mean_P=reshape(mean_P0,1,[])';
mean_P=mean_P(jk);

%% curve fitting
nihe_y=log(eta);
nihe_x=log(mean_P);
p=polyfit(nihe_x,nihe_y,1);
w1=exp(p(2));w2=p(1);
xi=linspace(min(mean_P),max(mean_P),20);
semilogy(xi,w1.*xi.^(w2),'linewidth',1.5,'color','k');

%% set gca
xlabel('Pressure (Pa)')
ylabel('\eta (Pa^-^1)')
grid on;
box on;
set(gca,'fontsize',16);

%% output
output=[mean_P,eta];
output_nihe=[w1,w2,min(mean_P),max(mean_P)];
save ./Out_sensitivity/Nur_1969_out.txt -ascii output
save ./Out_sensitivity/Nur_1969_nihe.txt -ascii output_nihe

%%
function [mean_P0,eta0]=process(P0,V0)
    jkf=find(~isnan(P0));
    P=P0(jkf);V=V0(jkf);
    figure(1)
    plot(P,V,'*-');hold on;
    eta0=zeros(1,length(P)-1);
    mean_P0=zeros(1,length(P)-1);
    for j=1:length(P)-1
        dV=V(j+1)-V(j);
        dP=P(j+1)-P(j);
        eta0(j)=dV/V(j)/dP;
        mean_P0(j)=(P(j+1)+P(j))/2;
    end
end

