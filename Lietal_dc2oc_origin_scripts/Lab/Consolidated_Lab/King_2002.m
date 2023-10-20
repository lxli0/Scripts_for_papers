clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=[2
5
10
15
20
30
40
50
60
90];
Velocity=[3.723412267	3.596360426	3.526487468	2.376285128	2.400999172	2.450427259	2.499842977	3.619266055	1.845183486	2.41146789	1.866513761
3.844334192	3.723644537	3.666485852	2.388753475	2.434074924	2.487609649	2.532924914	3.707568807	2.503440367	2.485779817	2.064678899
4.028968864	3.895575035	3.863824698	2.467249631	2.508408781	2.574325267	2.615521524	3.884174312	2.760321101	2.572477064	2.196788991
4.175511211	4.048469469	4.010336749	2.541577304	2.586880199	2.648677678	2.693992943	3.996559633	3.017201835	2.64266055	2.324770642
4.245798216	4.150506811	4.112404387	2.603560323	2.644756581	2.710660697	2.755975962	4.116972477	3.217889908	2.696330275	2.41146789
4.424494848	4.329203443	4.291101019	2.702818503	2.748133768	2.80579987	2.838745744	4.301605505	3.466743119	2.803669725	2.551834862
4.495195901	4.425333041	4.380868431	2.777381194	2.822684089	2.876237369	2.905064235	4.422018349	3.691513761	2.882110092	2.659174312
4.578641522	4.515100453	4.457931669	2.831330294	2.880764566	2.930198838	2.950781505	4.510321101	3.836009174	2.939908257	2.741743119
4.636658597	4.585831802	4.535005006	2.864702911	2.922362829	2.975916108	2.992392137	4.59059633	3.964449541	2.97293578	2.791284404
4.740816668	4.75354104	4.677315994	2.952445188	2.997754268	3.047169986	3.071896399	4.72706422	4.205275229	3.055504587	2.911009174
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
save ./Out_sensitivity/King_2002_out.txt -ascii output
save ./Out_sensitivity/King_2002_nihe.txt -ascii output_nihe

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

