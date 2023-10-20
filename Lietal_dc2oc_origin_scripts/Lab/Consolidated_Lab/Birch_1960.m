clc,clear
load('lajolla.mat')
lajolla=lajolla(20:150,:);
%% load data
Pressure=[10	500	1000	2.00E+03	4.00E+03	6.00E+03	1.00E+04];
Velocity=[5.6	nan	5.67	5.73	5.8	5.87	6
4.7	6.33	6.46	6.59	6.7	6.75	6.82
4.1	5.63	5.84	5.97	6.1	6.16	6.23
5.1	6.04	6.11	6.2	6.3	6.37	6.45
5	5.96	6.18	6.29	6.39	6.43	6.51
3.7	5.42	5.94	6.16	6.27	6.33	6.4
4.2	5.64	5.91	6.09	6.22	6.28	6.35
3.4	5.67	5.91	6.06	6.18	6.27	6.31
5.1	nan	5.95	6.07	6.22	6.28	6.37
5.6	nan	6.11	6.15	6.22	6.26	6.35
5.4	6.26	6.31	6.38	6.44	6.49	6.61
5.1	5.86	6.06	6.15	6.25	6.32	6.39
3.9	5	5.27	5.44	5.63	5.75	5.85
3.5	nan	4.73	5.02	5.38	5.58	5.89
5.9	nan	6.24	6.28	6.34	6.38	6.45
5.7	6.21	6.29	6.35	6.42	6.46	6.51
5.1	6.06	6.13	6.23	6.33	6.37	6.5
5.7	nan	6.42	6.46	6.51	6.55	6.61
6.1	6.28	6.33	6.37	6.43	6.48	6.57
5.4	5.63	5.76	5.87	5.98	6.04	6.13
5.7	6.13	6.19	6.25	6.3	6.34	6.41
6.4	nan	6.62	6.65	6.68	6.72	6.76
4.4	nan	6.27	6.35	6.43	6.48	6.56
5.4	nan	5.92	6.04	6.14	6.2	6.28
5.8	nan	6.02	6.08	6.15	6.21	6.31
5.49	nan	5.79	5.91	6.02	6.1	6.22
6.15	nan	6.24	6.3	6.36	6.4	6.46
4.4	nan	5.95	6.07	6.16	6.21	6.3
5.1	nan	6.33	6.43	6.49	6.54	6.6
6.73	nan	6.86	6.9	6.94	6.97	7.02
6.5	nan	6.97	7.01	7.05	7.07	7.1
5.7	nan	6.58	6.63	6.7	6.73	6.79
5.7	nan	6.43	6.48	6.53	6.57	6.64
6.4	nan	6.51	6.57	6.67	6.74	6.84
5.1	nan	6.43	6.52	6.6	6.64	6.71
5.7	6.92	6.98	7.05	7.13	7.16	7.21
4.8	nan	6.75	6.82	6.92	6.98	7.07
5.5	nan	6.46	6.53	6.6	6.65	6.71
4.9	nan	6.3	6.5	6.71	6.82	6.97
6.8	7.04	7.07	7.09	7.13	7.16	7.21
6.55	nan	6.64	6.67	6.71	6.75	6.82
6.14	nan	6.7	6.76	6.82	6.86	6.93
6.25	6.4	6.43	6.47	6.52	6.56	6.63
6.6	7.02	7.07	7.11	7.16	7.2	7.28
6	nan	6.37	6.46	6.55	6.64	6.79
6.4	6.67	6.72	6.76	6.81	6.84	6.91
6.76	nan	6.77	6.8	6.84	6.88	6.92
5.8	6.74	6.93	7.02	7.11	7.17	7.23
6.89	nan	7.17	7.21	7.27	7.31	7.35
7.6	nan	8.21	8.22	8.23	8.24	8.28
6.61	nan	7.2	7.32	7.41	7.47	7.54
7	nan	7.54	7.59	7.65	7.69	7.78
6.8	nan	7.73	7.79	7.88	7.93	8.01
7.5	7.69	7.75	7.8	7.86	7.92	8
7	7.82	7.89	8.01	8.13	8.19	8.28
7.42	nan	7.62	7.65	7.72	7.75	7.83
nan	7.4	7.49	7.6	7.75	7.85	8.02
7.7	nan	7.99	8.05	8.14	8.2	8.28
7.7	8.11	8.19	8.27	8.32	8.35	8.42
6.64	7.3	7.38	7.46	7.57	7.62	7.71
8.45	nan	8.67	8.69	8.72	8.75	8.78
6.6	7.49	7.56	7.65	7.79	7.85	7.92
6.9	7.74	7.78	7.81	7.85	7.9	7.95
7.17	7.65	7.68	7.73	7.79	7.82	7.87
5.2	nan	7.13	7.3	7.46	7.54	7.69
7.31	nan	7.69	7.81	7.89	7.94	8.01
6.3	nan	8.41	8.55	8.72	8.83	8.99
6.7	7.13	7.16	7.21	7.27	7.3	7.36
5.9	nan	7.81	7.91	7.99	8.01	8.07
];
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
    P0=1e5*Pressure;
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
save ./Out_sensitivity/Birch_1960_out.txt -ascii output
save ./Out_sensitivity/Birch_1960_nihe.txt -ascii output_nihe

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

