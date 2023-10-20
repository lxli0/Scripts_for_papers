clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=[0.1
5
10
15
20
30
40
50];
Velocity=[4.392	5.199	5.32	2.485	2.527	2.87
nan	5.754	6.008	nan	3.13	3.229
5.637	5.893	6.043	3.015	3.206	3.239
5.658	5.99	6.097	3.07	3.22	3.27
5.777	6.015	6.115	3.115	3.227	3.28
5.867	6.091	6.151	3.154	3.249	3.291
5.959	6.168	6.151	3.187	3.264	3.291
6.007	6.194	6.188	3.214	3.271	3.296];
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
save ./Out_sensitivity/Asef_2013_out.txt -ascii output
save ./Out_sensitivity/Asef_2013_nihe.txt -ascii output_nihe
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
