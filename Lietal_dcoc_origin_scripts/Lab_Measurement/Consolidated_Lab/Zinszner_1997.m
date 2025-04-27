clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Sample=[0.1000	1335.4681		0.1000	1748.4042		1.0000	2093.6939		0.1000	2515.2204		0.1000	3075.5651		0.1000	3307.7870
0.1857	1704.8153		0.4800	1879.1007		2.0465	2857.6071		0.4783	2590.6177		0.4740	3088.3039		0.4582	3324.2124
0.4705	1927.0057		1.0531	2020.0315		4.8777	3492.6355		1.9557	2754.3328		0.9940	3114.9007		1.0126	3343.4919
0.9662	2204.2714		1.9610	2205.0270		10.2843	3733.0631		4.8911	2932.2698		1.9521	3134.0543		1.9373	3347.8711
1.9573	2592.1217		4.9660	2508.3229		19.9276	3825.9350		9.9853	3069.4371		4.7276	3182.9247		4.4239	3356.1258
4.4817	3381.9461		10.3337	2759.7963		nan	nan		20.2528	3199.2242		9.7208	3198.4407		10.3702	3375.4683
10.3499	3773.6229		18.9887	2977.9573		nan	nan		45.9043	3395.4893		20.6538	3210.3050		nan	nan
20.4439	3954.9946		43.6125	3148.4300		nan	nan		nan	nan		nan	nan		nan	nan
40.9360	4033.1551		nan	nan		nan	nan		nan	nan		nan	nan		nan	nan
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
save ./Out_sensitivity/Zinszner_1997_out.txt -ascii output
save ./Out_sensitivity/Zinszner_1997_nihe.txt -ascii output_nihe
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

