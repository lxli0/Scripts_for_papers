clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Sample=[113.7085137	472.4997981	114.5716816	496.8491873	106.482626	478.1884348
152.0382395	523.0713115	189.8225826	571.7892353	152.612733	536.6152511
194.8773449	571.7044647	240.7989995	627.9942712	215.6897803	599.2154113
251.2445887	626.1829852	371.8812143	742.4860112	303.0182834	682.6822917
332.4134199	698.1799404	529.6653617	854.8960831	414.5284536	761.975828
463.1854257	795.496725	685.0220607	938.1628031	564.7769802	851.7027244
593.9574315	863.6623902	801.5395849	994.3678391	697.9913309	905.9561966
704.4372294	920.1448924	908.3473154	1038.082867	833.6308272	962.2963408
864.520202	998.0604082	1046.711875	1104.696243	915.9810937	995.6830929
1002.056277	1052.629789	1155.947054	1138.002931	1066.119122	1045.763221
1175.667388	1113.069777	1296.739063	1190.044631	1172.680134	1085.409989
1328.986291	1163.77001	1430.248726	1233.759659	1269.528933	1112.536725
1511.616162	1220.333277	1546.76625	1271.229683	1402.696758	1150.096822
1714.53824	1273.032443	1677.848465	1314.944711	1533.4627	1193.916934
1899.422799	1317.937787	1842.914958	1352.414735	1647.235177	1216.870326
2036.958874	1356.959904	1986.134415	1396.129763	1756.198072	1250.257078
2237.626263	1405.769731	2146.34601	1429.436451	1874.820839	1277.383814
2427.020202	1448.736714	2313.839952	1473.151479	2000.701596	1304.51055
2620.923521	1491.708745	2529.882861	1518.948175	2133.840343	1331.637286
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
    P0=1e3*Sample(:,2*i-1);
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
save ./Out_sensitivity/Linneman_2021_out.txt -ascii output
save ./Out_sensitivity/Linneman_2021_nihe.txt -ascii output_nihe
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

