clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Sample=[0.0000	3.3799	0.0169	1.8169
0.0089	3.5074	0.0404	2.0002
0.0260	3.6143	0.0489	2.0526
0.0355	3.6511	0.0962	2.2246
0.0571	3.7440	0.1969	2.3285
0.0980	3.9428	0.2947	2.3732
0.2075	4.1092	0.4054	2.4082
0.2514	4.1374	0.5144	2.4214
0.3139	4.1839	0.6052	2.4390
0.4062	4.2369	0.7143	2.4505
0.5167	4.2705	0.7997	2.4567
0.6202	4.2814	0.9052	2.4691
0.7169	4.3074	1.0087	2.4832
0.8227	4.3238	1.6035	2.5101
0.9192	4.3357	nan	nan
1.0203	4.3466	nan	nan
1.6247	4.3981	nan	nan
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
    P0=1e8*Sample(:,2*i-1);
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
save ./Out_sensitivity/Christensen_1985_out.txt -ascii output
save ./Out_sensitivity/Christensen_1985_nihe.txt -ascii output_nihe
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

