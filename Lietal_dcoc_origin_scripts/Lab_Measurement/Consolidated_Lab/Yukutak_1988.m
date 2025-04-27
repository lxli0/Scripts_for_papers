clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Sample=[0.0000	4.2359	0.0000	2.6272
1.6796	4.3218	2.9533	2.6571
8.6206	4.3971	6.7716	2.7224
9.6715	4.4722	10.5329	2.7345
19.7937	4.6657	20.8091	2.8300
27.2644	4.7625	27.8279	2.8837
45.4087	4.9938	46.6281	2.9384
59.8462	5.0694	61.1377	2.9868
80.1496	5.2365	82.0736	3.0357
100.4946	5.2480	100.8484	3.0667
120.2940	5.2971	120.1524	3.0918
153.4601	5.4217	153.3925	3.1297
199.4947	5.4825	201.6335	3.1747
251.9492	5.5651	253.0812	3.2140
300.1285	5.6153	302.3362	3.2058
349.9171	5.6548	351.0936	3.2331
402.3773	5.7159	404.6517	3.2430
];
Sample_number=size(Sample,2)/2;

%% create colormap
range=0:size(lajolla,1)-1;
range=range*(Sample_number-1)/(size(lajolla,1)-1)+1;
color=nan(Sample_number,3);
for i=1:3
    color(:,i)=interp1(range,lajolla(:,i),1:Sample_number);
end
%% Process
for i=1:Sample_number
    P0=1e6*Sample(:,2*i-1);
    V0=Sample(:,2*i);
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
save ./Out_sensitivity/Yukutak_1988_out.txt -ascii output
save ./Out_sensitivity/Yukutak_1988_nihe.txt -ascii output_nihe
%}
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

