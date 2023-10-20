clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=[20
50
80
100
150
200
250
300
400
500
600];
Velocity=[4.24	4.25	4.37	4.38	4.47	4.52	4.41	4.57	4.63	4.61	4.35	4.5	4.09	4.23	4.02	4.26	3.97	3.92	3.76	2.89	3.79	3.59	3.27	2.84	3.65
4.44	4.42	4.48	4.54	4.59	4.63	4.52	4.7	4.74	4.7	4.46	4.62	4.21	4.35	4.22	4.34	4.12	4.03	3.97	3.13	3.86	3.71	3.41	2.85	3.77
4.55	4.53	4.56	4.64	4.67	4.71	4.6	4.77	4.8	4.77	4.55	4.7	4.29	4.43	4.31	4.39	4.21	4.09	4.09	3.27	3.91	3.77	3.49	2.86	3.83
4.59	4.59	4.59	4.67	4.71	4.74	4.64	4.79	4.82	4.8	4.6	4.73	4.32	4.46	4.35	4.41	4.24	4.12	4.14	3.33	3.94	3.8	3.52	2.87	3.85
4.64	4.66	4.64	4.72	4.75	4.78	4.69	4.82	4.86	4.83	4.68	4.78	4.37	4.5	4.4	4.45	4.27	4.14	4.2	3.43	3.97	3.84	3.56	2.88	3.88
4.67	4.71	4.66	4.75	4.78	4.8	4.71	4.84	4.89	4.85	4.72	4.8	4.4	4.52	4.42	4.46	4.3	4.16	4.24	3.48	3.98	3.85	3.58	2.88	3.9
4.68	4.73	4.68	4.77	4.79	4.82	4.73	4.85	4.91	4.86	4.75	4.82	4.43	4.53	4.44	4.48	4.31	4.17	4.27	3.52	3.99	3.86	3.59	2.89	3.91
4.7	4.74	4.7	4.78	4.8	4.83	4.75	4.86	4.92	4.88	4.77	4.84	4.45	4.54	4.45	4.49	4.33	4.18	4.29	3.55	4	3.87	3.61	2.89	3.92
4.72	4.77	4.72	4.81	4.83	4.85	4.77	4.88	4.94	4.89	4.8	4.86	4.48	4.56	4.46	4.51	4.35	4.2	4.33	3.59	4.02	3.89	3.63	2.9	3.94
4.74	4.78	4.74	4.83	4.84	4.87	4.79	4.89	4.95	4.9	4.83	4.88	4.5	4.57	4.48	4.53	4.37	4.21	4.36	3.62	4.03	3.9	3.64	2.9	3.95
4.76	4.79	4.75	4.84	4.85	4.88	4.8	4.9	4.96	4.92	4.84	4.89	4.52	4.58	4.49	4.54	4.39	4.23	4.38	3.64	4.04	3.91	3.65	2.91	3.96
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
save ./Out_sensitivity/Wang_2005_out.txt -ascii output
save ./Out_sensitivity/Wang_2005_nihe.txt -ascii output_nihe
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
