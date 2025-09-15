clc,clear
close all
%{
% Create figure with normalized units and specified position
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.3, 0.2]);

% Scatter plot of data Dm
scatter(1:length(Dm), Dm+Mc, 30, 'k', 'filled');

% Labeling the axes
xlabel('Event Index', 'FontSize', 16);
ylabel('Magnitude', 'FontSize', 16);

% Formatting the axes
set(gca, 'FontSize', 16);
grid on;
grid minor;
box on;

% Title
title('Otaniemi', 'FontSize', 18, 'FontWeight', 'bold');

%}

%% Load data
D=load('Otaniemi_Mw.txt'); % only Decatur: (2,3):(northing,easting) and Otaniemi: (2,3):(easting,northing); other lat,lon
time_flag=1;
timeleg=['Days Since 2018/6/4'];
figure('units','normalized','position',[0.1,0.1,0.7,0.9])
MT_xlim=[0 50];
MT_ylim=[-0.5 2];
b_ylim=[0.8 1.6];
m_croner_ylim=[1.2 2];

D=sortrows(D,1);
T0=D(:,1);
if time_flag==2
    T0=1:1:length(T0);
end
X=D(:,2);
Y=D(:,3);
X=X-mean(X);
Y=Y-mean(Y);
M0=D(:,5);

diffM0=diff(M0);
nonzero_elements=diffM0(diffM0~=0);
delta= min(abs(nonzero_elements))/2;

%% map
subplot(3,3,1)
scatter(X,Y,30,T0,'filled');
colormap(summer);
cb=colorbar; 
if time_flag==2
    title(cb, '#');
else
    title(cb,timeleg);
end
xlabel('X (m)');
ylabel('Y (m)');
grid on;box on;grid minor;set(gca,'fontsize',16)
axis equal


%% MFD and DMFD
dm=0.1;
% get confidenc einterval of parameters
boostrap=1000;
seeds = 1:boostrap; 
for k=1:boostrap
    rng(seeds(k));
    random_Magn = M0(randi(length(M0),1,length(M0)));
    m=floor(min(random_Magn)/dm)*dm:dm:max(random_Magn);
    n0=hist(random_Magn,m);   
    [~,ind]=max(n0);
    Mc_MLE=m(ind)+0.2;
    [Mc_MBS,b_es,number]=MBS_MLE_discrete(random_Magn,delta);
    Mc=max(Mc_MLE,Mc_MBS);
    jkf=find(random_Magn-Mc>=-1e-11);
    Dm=random_Magn(jkf);
    Dm=Dm-Mc;
    Mmin=0;
    [b_GR_bt(k),loglikelihood_GR]=Estimation_GR_discrete(Dm,Mmin,delta);
    [b_TGR_bt(k),mcorner_TGR_bt(k),loglikelihood_TGR]=Estimation_TGR_gridsearch_discrete(Dm,Mmin,delta,b_GR_bt(k));
    mcorner_TGR_bt(k)=Mc+mcorner_TGR_bt(k);
    R=loglikelihood_TGR-loglikelihood_GR;
    P_LR_bt(k)=1-chi2cdf(2*R,1);
end
b_GR_percentiles = prctile(b_GR_bt, [5, 95])
b_TGR_percentiles = prctile(b_TGR_bt, [5, 95])
mcorner_TGR_percentiles = prctile(mcorner_TGR_bt, [5, 95])
%}
% MLE 
m=floor(min(M0)/dm)*dm:dm:max(M0);
n0=hist(M0,m);   
[~,ind]=max(n0);
Mc_MLE=m(ind)+0.2;
[Mc_MBS,b_es,number]=MBS_MLE_discrete(M0,delta);
Mc=max(Mc_MLE,Mc_MBS);
jkf=find(M0-Mc>=-1e-11);
Dm=M0(jkf);
Dt=T0(jkf);
Dm=Dm-Mc;
Mmin=0;
[b_GR,loglikelihood_GR]=Estimation_GR_discrete(Dm,Mmin,delta);
[b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR_gridsearch_discrete(Dm,Mmin,delta,b_GR);
R=loglikelihood_TGR-loglikelihood_GR;
P_LR=1-chi2cdf(2*R,1);

subplot(3,3,2)
plt_single_MFD(M0,Dm,Mc,b_GR,b_TGR,mcorner_TGR,P_LR);

subplot(3,3,3)

%% EQ rate
color_ope=[0 0 0.5;0.4 0.4 0.4];
[hh,Vtotalmax]=plt_Otaniemi_operation(color_ope);

color_mc=[ 0.8 0.2 0;1 0.6 0.2];
t4count=0:0.1:50;
for i=1:length(t4count)
    jkf=find(T0<=t4count(i));
    M_count=M0(jkf);
    EQ_count(i,1)=length(M_count);
    EQ_count(i,2)=max([length(M_count(M_count>=Mc-1e-6)),0]);
end

yyaxis right
for j=1:2
    hh(2+j)=plot(t4count,EQ_count(:,j)/max(EQ_count(:,j)),'-','color',color_mc(j,:),'linewidth',2.5);hold on;
end
addstage
ylabel('Normalized Cumulative Values')
% Create legend labels
label1 = 'Pressure';
label2 = ['Injected Volume / ', num2str(Vtotalmax, '%.0f'), ' m³']; % Space around slash
label3 = ['EQ Count / ', num2str(max(EQ_count(:,1)), '%.0f')]; % Fixed parenthesis
label4 = ['EQ Count (m ≥ 0.4) / ', num2str(max(EQ_count(:,2)), '%.0f')]; % Rephrased for clarity
ylim([0 1])
% Create legend
legend(hh, label1, label2, label3, label4, ...
    'NumColumns', 2, ...
    'Location', 'Northwest');
% Adjust plot appearance
xlim([0 50]);
grid on; box on; grid minor;
set(gca, 'FontSize', 16); % Increase font size for readability


figure('units','normalized','position',[0.1,0.1,0.7,0.9])
%% temporal evolution
Ntotal=length(Dm);
min_size=100;
max_size=1e3;
w0=ceil(Ntotal/max_size)+1:1:min([ceil(Ntotal/min_size)-1,ceil(Ntotal/max_size)+1+4]); % number of windows
[P_LR_median_eachwindow,bGR_median_eachnode,bTGR_median_eachnode,mcorner_median_eachnode]=get_temporal_MFD(Dm,delta,Mmin,min_size,max_size,w0);
mcorner_median_eachnode=mcorner_median_eachnode+Mc;
DmplusMc=Dm+Mc;

subplot(3,2,3)
median_PLR=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,P_LR_median_eachwindow,1,time_flag,timeleg);
xlim([MT_xlim(1),MT_xlim(2)]);
yyaxis left
ylim([MT_ylim(1),MT_ylim(2)]);

subplot(3,2,4)
median_bGR=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,bGR_median_eachnode,2,time_flag,timeleg);
xlim([MT_xlim(1),MT_xlim(2)]);
yyaxis left
ylim([MT_ylim(1),MT_ylim(2)]);
yyaxis right
ylim([b_ylim(1),b_ylim(2)]);
addstage

subplot(3,2,5)
median_bTGR=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,bTGR_median_eachnode,3,time_flag,timeleg);
xlim([MT_xlim(1),MT_xlim(2)]);
yyaxis left
ylim([MT_ylim(1),MT_ylim(2)]);
yyaxis right
ylim([b_ylim(1),b_ylim(2)]);
addstage

subplot(3,2,6)
median_mcorner=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,mcorner_median_eachnode,4,time_flag,timeleg);
xlim([MT_xlim(1),MT_xlim(2)]);
yyaxis left
ylim([MT_ylim(1),MT_ylim(2)]);
yyaxis right
ylim([m_croner_ylim(1),m_croner_ylim(2)]);
addstage


function [P_LR_median_eachwindow,bGR_median_eachnode,bTGR_median_eachnode,mcorner_median_eachnode]=get_temporal_MFD(x,delta,mmin,min_size,max_size,w0)
    Ntotal=length(x);
    num_rand=100; % for each w0, number of random partition
    P_LR_median_eachwindow=zeros(length(w0),Ntotal);
    bGR_median_eachnode=zeros(length(w0),Ntotal);
    bTGR_median_eachnode=zeros(length(w0),Ntotal);
    mcorner_median_eachnode=zeros(length(w0),Ntotal);
    
    for k=1:length(w0)
        w=w0(k);
        P_LR=zeros(num_rand,Ntotal);
        bGR=zeros(num_rand,Ntotal);
        bTGR=zeros(num_rand,Ntotal);
        mcorner=zeros(num_rand,w+1);
        node=zeros(num_rand,w+1);
        for i=1:num_rand
            node(i,:)=generateNumbers(w,length(x),min_size,max_size);% w+1 nodes, including the first and last events
            for j=1:w
                xx=x(node(i,j):node(i,j+1));
                [b_GR,loglikelihood_GR]=Estimation_GR_discrete(xx,mmin,delta);
                [b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR_gridsearch_discrete(xx,mmin,delta,b_GR);
                R=loglikelihood_TGR-loglikelihood_GR;
                P_LR(i,node(i,j):node(i,j+1))=1-chi2cdf(2*R,1);
                bGR(i,node(i,j):node(i,j+1))=b_GR;    
                bTGR(i,node(i,j):node(i,j+1))=b_TGR;       
                mcorner(i,node(i,j):node(i,j+1))=mcorner_TGR;       
            end
        end
        P_LR_median_eachwindow(k,:)=median(P_LR,1); 
        bGR_median_eachnode(k,:)=median(bGR,1); 
        bTGR_median_eachnode(k,:)=median(bTGR,1); 
        mcorner_median_eachnode(k,:)=median(mcorner,1);         
    end
end

function nonde_order=generateNumbers(n,Ntotal,min_size,max_size) % n is the window number
    nonde_order=zeros(1,n+1);
    target_sum=Ntotal-n*min_size; % available number for arrangement
    scaled_numbers=1e5*ones(1,n);
    while max(scaled_numbers)>=max_size-min_size
        random_numbers=rand(1,n);
        scaled_numbers=floor(random_numbers/sum(random_numbers)*target_sum); % number of events (-100) for each window
    end
    nonde_order(1)=1;
    for i=2:n
        nonde_order(i)=nonde_order(i-1)+min_size+scaled_numbers(i-1);
    end
    nonde_order(n+1)=Ntotal;
end

function [meadian_plt_variable]=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,plt_variable,plt_variable_flag,time_flag,timeleg)
    yyaxis left
    scatter(T0,M0,30,[0.75 0.75 0.75],'filled');
    hold on;
    scatter(Dt,DmplusMc,30,'k','filled');
    ylabel('Magntidue');
    if time_flag==1
        xlabel(timeleg)
    elseif time_flag==2
        xlabel('#')
    end
    ax = gca; 
    ax.YColor = 'k';

    yyaxis right
    col=autumn(length(w0));
    for k=1:length(w0)
        h(k)=plot(Dt,plt_variable(k,:),'-','color',col(k,:),'linewidth',2.5);hold on;
        leg{k}=['$w=$',num2str(w0(k))];
    end
    meadian_plt_variable=median(plt_variable,1);
    h(k+1)=plot(Dt,meadian_plt_variable,'-b','linewidth',2.5);hold on;
    
    leg{k+1}=['median'];
    if plt_variable_flag==1
        ylabel('$p$','Interpreter','latex')
        ylim([0 1]);
        addstage
        legend(h,leg, 'Interpreter', 'latex','location','northwest','Numcolumns',3);
    elseif plt_variable_flag==2
        ylabel('$b_{GR}$','Interpreter','latex')
    elseif plt_variable_flag==3
        ylabel('$b_{TGR}$','Interpreter','latex')
    elseif plt_variable_flag==4
        ylabel('$m_{corner}$','Interpreter','latex')
    end  
    ax = gca; 
    ax.YColor = 'k';
    set(gca,'fontsize',16);
    grid on;box on;grid minor;
end

function []=addstage()
    % Shade Stages
    colors = {'r', 'g', 'm', 'c', 'y'};
    t_stages = [0 13 23 27 36.5 50.5];
    alpha_plot = 0.1;
    for i = 1:length(t_stages) - 1
        x = [t_stages(i), t_stages(i+1), t_stages(i+1), t_stages(i)];
        y = [0, 0, 100, 100];
        fill(x, y, colors{i}, 'FaceAlpha', alpha_plot, 'EdgeColor', 'none');
    end
end
%}