clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=[0	0.1	0.4	0.7	1	2	3	5	10];
Velocity=[3.64	4.336	5.181	5.524	5.774	6.029	6.138	6.235	6.368
2.421	3.008	3.287	3.447	3.535	3.64	3.705	3.733	3.758
5.418	5.631	5.806	5.903	5.959	6.014	6.071	6.133	6.236
3.182	3.482	3.803	3.886	3.92	3.96	3.977	3.987	3.994
4.15	4.824	5.383	5.569	5.706	5.827	5.954	6.064	6.208
2.686	3.016	3.251	3.332	3.382	3.476	3.531	3.566	3.593
6.299	6.35	6.445	6.48	6.496	6.529	6.562	6.596	6.599
3.366	3.373	3.378	3.381	3.384	3.388	3.389	3.389	3.386
6.331	6.334	6.344	6.353	6.362	6.391	6.414	6.438	6.485
3.095	3.105	3.107	3.112	3.117	3.131	3.143	3.168	3.222
6.65	6.666	6.72	6.756	6.783	6.834	6.886	6.922	7.018
3.724	3.753	3.764	3.774	3.781	3.795	3.907	3.823	3.843
5.876	5.892	5.911	5.921	5.927	5.942	5.957	5.984	6.007
3.758	3.756	3.75	3.743	3.737	3.717	3.697	3.648	3.528
];
Sample_number=size(Velocity,1);

%% create colormap
range=0:size(lajolla,1)-1;
range=range*(Sample_number-1)/(size(lajolla,1)-1)+1;
color=nan(Sample_number,3);
for i=1:3
    color(:,i)=interp1(range,lajolla(:,i),1:Sample_number);
end
%% Process
for i=1:Sample_number
    P0=1e8*Pressure;
    V0=Velocity(i,:);
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
save ./Out_sensitivity/Simmons_1965_out.txt -ascii output
save ./Out_sensitivity/Simmons_1965_nihe.txt -ascii output_nihe

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

