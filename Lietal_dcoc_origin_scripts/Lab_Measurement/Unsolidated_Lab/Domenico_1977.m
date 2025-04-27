clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=[0.3
0.4
0.5
0.77
1
1.5
2
3
4
4.5
5];
Velocity=0.3048*[nan	nan	nan	nan	nan	nan	1.801595738	nan
6.233287707	2.844679246	1.864814786	nan	6.431007645	2.748610897	1.926854734	nan
6.190504392	3.014309041	2.004606842	nan	6.480988773	2.9553937	nan	nan
6.402133087	3.310981421	2.217801643	2.180895995	6.549421961	3.274834756	2.133066013	nan
6.486392357	3.522740851	2.416288658	2.290841233	6.617895971	3.55656884	2.258025638	2.069817364
6.580489501	3.935441319	2.658127326	2.584361426	6.73599731	3.944443084	2.539229022	2.319641362
6.780744175	4.220739677	2.855729151	2.745057603	6.854098649	4.162874088	2.682422983	2.544440776
6.883927463	4.653148577	3.177280386	3.029725888	7.021351046	4.562246849	2.93748595	2.736704149
7.093366331	5.043199055	3.395604754	3.240650969	7.106957166	4.848634771	3.136076909	2.922776013
7.187365424	5.24338836	3.452937882	nan	7.231386091	5.004375043	3.229058411	3.015771121
7.249595699	5.443741085	3.503121392	3.31133005	7.343214275	5.135063516	3.296988112	3.114998562
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
    P0=6895*1e3*Pressure;
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
save ./Out_sensitivity/Domenico_1977_out.txt -ascii output
save ./Out_sensitivity/Domenico_1977_nihe.txt -ascii output_nihe

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

