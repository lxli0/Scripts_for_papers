clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=[250.0000 
500.0000 
750.0000 
1000.0000 
1500.0000 
2000.0000 
3000.0000 
4000.0000 
5000.0000 ];
Velocity=[3.3175 	3.7461 	3.0708 	3.4149 	3.3890 	2.9929 	3.6097 
3.4727 	3.9532 	3.1806 	3.5701 	3.3818 	3.1416 	3.7130 
3.5695 	4.0435 	3.3293 	3.6733 	3.4591 	3.2643 	3.8033 
3.6533 	4.0883 	3.3676 	3.7442 	3.5429 	3.3480 	3.9000 
3.7623 	4.1909 	3.5676 	3.8208 	3.6974 	3.4896 	3.9896 
3.8844 	4.2546 	3.6506 	3.8649 	3.7935 	3.6442 	4.0986 
4.0506 	4.3623 	3.9013 	3.9403 	3.9598 	3.8429 	4.3169 
4.2818 	4.4507 	4.1065 	4.0675 	4.1519 	4.0415 	4.5091 
4.4351 	4.5974 	4.2338 	4.1559 	4.2663 	4.1753 	4.6623  
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
    P0=6895*Pressure;
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
save ./Out_sensitivity/Smith_2010_out.txt -ascii output
save ./Out_sensitivity/Smith_2010_nihe.txt -ascii output_nihe

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

