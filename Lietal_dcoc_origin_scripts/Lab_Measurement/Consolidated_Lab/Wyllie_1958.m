clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Sample=[0.1135	9.7466		0.1095	11.2206		0.1162	10.0952			0.0850	12.1453		0.0315	11.3265
0.6000	10.9835		0.8947	11.2964		1.0472	11.1256			1.7031	12.5277		0.9796	12.1424
1.0647	11.3947		1.7907	11.7762		2.1131	11.6915			2.5301	12.8148		2.0367	12.4198
1.5592	11.6837		2.5779	12.1088		3.0535	11.9452			3.2944	13.0033		2.9924	12.6347
2.5772	12.0232		3.4257	12.3316		4.1014	12.1686			4.1613	13.1714		3.9628	12.7045
3.5487	12.3198		4.2267	12.4625		5.0712	12.2388			5.0295	13.2358		4.9966	12.8108
4.5194	12.5062		5.1048	12.6059		6.1023	12.2848			5.8602	13.2328		5.9926	12.8704
5.5355	12.6072		5.8744	12.6572		7.0254	12.2937			nan	nan		6.9762	12.9041
6.5516	12.6898		nan	nan		8.0258	12.3396			nan	nan		7.9602	12.9066
7.5213	12.7601		nan	nan		9.0719	12.3367			nan	nan		nan	nan
8.5525	12.8121		nan	nan		10.0413	12.3519			nan	nan		nan	nan
10.0294	12.8105		nan	nan		nan	nan			nan	nan		nan	nan
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
    P0=1e3*6895*Sample(:,2*i-1);
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
save ./Out_sensitivity/Wyllie_1958_out.txt -ascii output
save ./Out_sensitivity/Wyllie_1958_nihe.txt -ascii output_nihe
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

