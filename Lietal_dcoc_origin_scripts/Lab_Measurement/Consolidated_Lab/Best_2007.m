clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=[5
10
20
40
60];
Velocity=[3266	nan	3842	2218	2001	3763	2198	4173	2440	2257	5981	3057	nan	nan	nan	3839	2212	3745	2088	2077	3427	2235	3499	3372	4063	3578	4145
3415	nan	nan	nan	nan	3851	2298	nan	nan	nan	5898	3052	5972	3112	3096	3996	2339	nan	nan	nan	3460	2205	3578	3379	nan	3611	nan
3554	2140	nan	nan	nan	3949	2395	nan	nan	nan	5890	3058	nan	nan	nan	4134	2469	nan	nan	nan	3498	2138	3652	3401	nan	3650	nan
3704	2230	4299	2596	2387	4044	2471	4374	2670	2563	5912	3071	5996	3123	3115	4236	2556	4203	2523	2492	3541	2083	3704	3448	4278	3690	4241
3796	2317	nan	nan	nan	4101	2519	nan	nan	nan	5941	3080	nan	nan	nan	4281	2584	nan	nan	nan	3509	2012	3671	3469	nan	3684	nan
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
save ./Out_sensitivity/Best_2007_out.txt -ascii output
save ./Out_sensitivity/Best_2007_nihe.txt -ascii output_nihe

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

