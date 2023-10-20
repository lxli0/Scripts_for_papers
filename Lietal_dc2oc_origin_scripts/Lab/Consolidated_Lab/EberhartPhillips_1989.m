clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Sample=[0.0210	4.7795	0.0205	2.9374		0.0091	3.8754	0.0098	2.1063		0.0413	3.7148	0.0432	2.1824		0.0213	2.7152	0.0706	1.5417		0.0284	4.0778	0.0297	2.3672		0.0178	3.1512	0.0237	1.7403
0.0288	4.7966	0.0305	2.9544		0.0504	4.0304	0.0512	2.2788		0.0637	3.7847	0.0634	2.2702		0.0418	2.8725	0.0876	1.6377		0.0485	4.2260	0.0497	2.4891		0.0289	3.2312	0.0403	1.8023
0.0422	4.7870	0.0428	2.9800		0.0983	4.1155	0.1014	2.3899		0.0815	3.8282	0.0835	2.3047		0.0510	2.9949	0.1000	1.6375		0.0585	4.3309	0.0587	2.5852		0.0411	3.2487	0.0569	1.8733
0.0499	4.8041	0.0506	2.9970		0.0973	4.1678	0.1883	2.4907		0.1017	3.8715	0.1036	2.3480		0.0623	3.0035	0.1393	1.7070		0.0685	4.3480	0.0698	2.6023		0.0533	3.3109	0.0691	1.8818
0.0611	4.8122	0.0594	2.9965		0.1897	4.2161	0.2005	2.4901		0.1530	3.9487	0.1515	2.4075		0.0691	3.0647	0.1528	1.7155		0.0908	4.4172	0.0909	2.6979		0.0611	3.2928	0.0913	1.9169
0.0722	4.8203	0.0717	3.0221		0.2020	4.2417	0.2874	2.5125		0.2043	3.9902	0.2017	2.4580		0.0894	3.0994	0.1877	1.7675		0.1019	4.4342	0.1032	2.6886		0.0699	3.3640	0.1001	1.9166
0.0888	4.8194	0.0906	3.0298		0.2877	4.2728	0.2963	2.5208		0.2912	4.0563	0.2897	2.4973		0.1030	3.1692	0.2012	1.7761		0.1397	4.5115	0.1388	2.7835		0.0888	3.3455	0.1378	1.9512
0.0988	4.8539	0.1006	3.0468		0.2989	4.2985	0.3909	2.5515		0.3024	4.0823	0.3020	2.4966		0.1412	3.2036	0.2922	1.8535		0.1520	4.5373	0.1510	2.7829		0.1032	3.3630	0.1500	1.9419
0.1388	4.8428	0.1406	3.0532		0.3891	4.3207	0.4009	2.5511		0.3947	4.0679	0.3933	2.5179		0.1525	3.2472	0.3045	1.8708		0.1898	4.5707	0.1921	2.8425		0.1386	3.3976	0.1910	1.9854
0.1489	4.8860	0.1495	3.0527		0.3968	4.3291	0.4911	2.5820		0.4070	4.0939	0.4022	2.5352		0.1917	3.2554	0.3909	1.8783		0.2009	4.5702	0.2010	2.8422		0.1475	3.4509	0.2020	1.9940
0.1877	4.8662	0.1895	3.0679		0.4904	4.3598	nan	nan		0.4939	4.1155	0.4935	2.5476		0.2053	3.2990	0.4033	1.8869		0.2909	4.6452	0.2888	2.8909		0.1906	3.4229	0.2895	2.0181
0.2022	4.9005	0.1984	3.0673		nan	nan	nan	nan		nan	nan	nan	nan		0.2940	3.3415	0.4941	1.9118		0.2998	4.6448	0.3010	2.8729		0.2029	3.4939	0.2995	2.0536
0.2866	4.8780	0.2884	3.0796		nan	nan	nan	nan		nan	nan	nan	nan		0.3030	3.3588	nan	nan		0.3887	4.6497	0.3911	2.8952		0.2870	3.5004	0.3903	2.0509
0.2989	4.8948	0.2950	3.0792		nan	nan	nan	nan		nan	nan	nan	nan		0.3894	3.3838	nan	nan		0.4021	4.6842	0.4044	2.9122		0.3003	3.5535	0.4013	2.0684
0.3878	4.9071	0.3873	3.1002		nan	nan	nan	nan		nan	nan	nan	nan		0.4051	3.3924	nan	nan		0.4898	4.6804	0.4933	2.9346		0.3878	3.5331	0.4888	2.0926
0.3978	4.9241	0.3984	3.1083		nan	nan	nan	nan		nan	nan	nan	nan		0.4948	3.4086	nan	nan		nan	nan	nan	nan		0.4011	3.5862	nan	nan
0.4866	4.9365	0.4895	3.1293		nan	nan	nan	nan		nan	nan	nan	nan		nan	nan	nan	nan		nan	nan	nan	nan		0.4907	3.5925	nan	nan
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
save ./Out_sensitivity/EberhartPhillips_1989_out.txt -ascii output
save ./Out_sensitivity/EberhartPhillips_1989_nihe.txt -ascii output_nihe
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

