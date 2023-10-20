clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Sample=[60	6.326	3.343		60	5.728	3.158		60	6.185	3.187		60	5.784	3.143		60	6.329	3.474		30	6.468	3.334		60	6.379	3.375		60	6.312	3.361		60	6.581	3.366		30	nan	3.333		60	6.522	3.417		60	6.607	nan		60	6.515	3.387		60	nan	3.373
30	6.22	3.317		30	5.19	2.979		30	6.08	3.166		30	5.325	2.998		30	6.193	3.337		15	6.405	3.239		15	6.049	3.275		30	6.175	3.331		30	6.533	3.33		5	nan	3.095		30	6.491	3.356		30	6.579	nan		30	6.391	3.325		30	nan	3.338
15	6.076	3.268		15	4.562	2.737		15	5.633	3.033		15	4.79	2.79		15	6.043	3.2875		7.5	6.293	3.145		5	5.599	3.119		15	6.014	3.281		15	6.435	3.278		nan	nan	nan		15	6.328	3.232		15	6.504	nan		15	6.193	3.25		14	nan	3.31
10	5.987	3.238		10	4.22	2.601		10	5.27	nan		10	4.508	2.672		10	5.943	3.246		5	6.252	3.079		nan	nan	nan		10	5.906	3.25		10	6.347	3.228		nan	nan	nan		10	6.186	3.145		10	6.446	nan		10	6.056	3.186		10	nan	3.228
5	5.824	3.176		5	3.72	2.37		5	4.43	nan		5	4.08	2.474		5	5.681	3.148		2.5	6.202	3.01		nan	nan	nan		5	5.685	3.188		5	6.11	3.122		nan	nan	nan		5	5.928	2.965		5	6.261	nan		5	5.765	3.027		5	nan	3.142
nan	nan	nan		nan	nan	nan		5.7	nan	2.76		nan	nan	nan		nan	nan	nan		nan	nan	nan		nan	nan	nan		nan	nan	nan		nan	nan	nan		nan	nan	nan		nan	nan	nan		nan	nan	nan		nan	nan	nan		nan	nan	nan];
Sample_number=size(Sample,2)/3*2;

%% create colormap
range=0:size(lajolla,1)-1;
range=range*(Sample_number-1)/(size(lajolla,1)-1)+1;
color=nan(Sample_number,3);
for i=1:3
    color(:,i)=interp1(range,lajolla(:,i),1:Sample_number);
end
%% Process
for i=1:Sample_number
    judge=rem(i,2)==0;
    if judge==0
        P0=1e6*Sample(:,3*(i+1)/2-2);
        V0=Sample(:,3*(i+1)/2-1);
    else
        P0=1e6*Sample(:,3*i/2-2);
        V0=Sample(:,3*i/2);
    end   
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
save ./Out_sensitivity/Peacock_1994_out.txt -ascii output
save ./Out_sensitivity/Peacock_1994_nihe.txt -ascii output_nihe
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

