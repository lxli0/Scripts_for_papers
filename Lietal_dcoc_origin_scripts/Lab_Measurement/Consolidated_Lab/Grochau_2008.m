clc,clear

load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=[3
7
10
14
17
21
28
34
41];
Velocity=[2362.391831	2625.891312	2456.282451	2598.632745	2592.575286	2892.419522
2622.862582	2940.879197	2744.011769	2886.362063	2880.304604	3007.51125
2831.844929	3022.654898	2919.678089	3034.769817	3028.712357	3216.493596
2865.160955	3110.488058	3001.45379	3152.890273	3140.775355	3361.87262
2910.5919	3237.694704	3083.229491	3271.01073	3231.637245	3446.677051
2977.223953	3274.03946	3143.804084	3331.585324	3307.355486	3504.222915
3031.741087	3355.815161	3231.637245	3434.562132	3404.274836	3604.170993
3074.143302	3401.246106	3271.01073	3492.107996	3486.050537	3664.745587
3113.516788	3443.648321	3325.527864	3540.56767	3528.452752	3701.090343
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
    P0=1e6*Pressure;
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
save ./Out_sensitivity/Grochau_2008_out.txt -ascii output
save ./Out_sensitivity/Grochau_2008_nihe.txt -ascii output_nihe

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

