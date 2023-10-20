clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=10^6*[1
2
5
10
20
30
50
70
90];
Velocity=[nan	nan	5415	5238	5464	5409	4092	nan	4441	4168	nan	nan	nan	nan	2494	nan	nan	nan	nan	nan	2287	nan
nan	5344	5446	5277	5512	5458	4138	3831	4533	4271	4216	nan	nan	nan	2524	nan	nan	nan	nan	nan	2318	nan
nan	5353	5471	5305	5549	5531	4257	4120	4589	4534	4427	nan	nan	2946	2575	2437	nan	nan	nan	nan	2352	2163
nan	5377	5493	5382	5554	5601	4487	4462	4627	4793	4557	nan	2935	2970	2629	2577	nan	nan	nan	nan	2413	2153
nan	5392	5504	5447	5597	5658	4577	4746	4718	4915	4645	nan	3004	2982	2667	2820	nan	nan	2559	nan	2461	2221
nan	5443	5555	5495	5642	5685	4774	4720	4780	4984	4910	nan	3053	3018	2696	3025	nan	nan	2595	nan	2499	2213
nan	5454	5587	5512	5628	5711	4825	4885	4787	5018	4962	nan	3061	3086	2752	3050	nan	nan	2759	nan	2621	2260
nan	5466	5598	5528	5646	5721	4784	4912	4832	5052	5204	nan	3058	3134	2843	3062	nan	nan	2729	nan	2658	2262
nan	5478	5606	5542	5637	5723	4846	5011	4971	5091	5050	nan	3061	3140	2896	3083	nan	nan	2849	nan	2709	2304
]';
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
    P0=Pressure;
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
save ./Out_sensitivity/Miller_2021_out.txt -ascii output
save ./Out_sensitivity/Miller_2021_nihe.txt -ascii output_nihe
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

