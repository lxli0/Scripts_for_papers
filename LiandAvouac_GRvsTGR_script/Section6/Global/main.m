clc,clear
close all
%% Load data
D=load('Global.txt'); % only Decatur: (2,3):(northing,easting) and Otaniemi: (2,3):(easting,northing); other lat,lon
jkf=find(D(:,1)>=1920);
D=D(jkf,:);
time_flag=1;
timeleg=['Year'];
figure('units','normalized','position',[0.1,0.1,0.7,0.9])
MT_xlim=[1920 2020];
MT_ylim=[5 10];
b_ylim=[0.95 1.2];
m_croner_ylim=[7.5 10];

D=sortrows(D,1);
T0=D(:,1);
if time_flag==2
    T0=1:1:length(T0);
end
Lat=D(:,2);
Lon=D(:,3);
M0=D(:,5);

% one can reshuffle to earthquakes to see check if temporal variaitons are significant
%randomIndex = randperm(length(M0));
%M0 = M0(randomIndex);

diffM0=diff(M0);
nonzero_elements=diffM0(diffM0~=0);
delta= min(abs(nonzero_elements))/2;

%% map
%{
scatter(Lon,Lat,30,T0,'filled');
colormap(summer);
cb=colorbar; 
if time_flag==2
    title(cb, '#');
else
    title(cb,timeleg);
end
xlabel('Longitude (°)');
ylabel('Latitude (°)');
grid on;box on;grid minor;set(gca,'fontsize',16)
axis equal
ylim([-90 90])
xlim([-180 180])
%}
%% MFD and DMFD
dm=0.1;
% get confidenc einterval of parameters
boostrap=1000;
for k=1:boostrap
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

subplot(3,3,1)
plt_single_MFD(M0,Dm,Mc,b_GR,b_TGR,mcorner_TGR,P_LR);

    
%{
subplot(3,3,3)
diffm=diff(M0);
diffm0=diffm(diffm>=-1e-11);
[Mc_prime,b_pos_1,number_abovemc]=MBS_MLE_discrete(diffm0,delta);
diffm=diffm0(diffm0>=Mc_prime-1e-11);
diffm=diffm-Mc_prime;
[b_pos_2,loglikelihood_GR_dm]=Estimation_GR_discrete(diffm,0,delta);
plt_single_MFD_DM(diffm0,diffm,Mc_prime,b_pos_2)
%}
%% temporal evolution
Ntotal=length(Dm);
min_size=100;
max_size=1e3;
w0=ceil(Ntotal/max_size)+1:1:min([ceil(Ntotal/min_size)-1,ceil(Ntotal/max_size)+1+4]); % number of windows
[P_LR_median_eachwindow,bGR_median_eachnode,bTGR_median_eachnode,mcorner_median_eachnode]=get_temporal_MFD(Dm,delta,Mmin,min_size,max_size,w0);
mcorner_median_eachnode=mcorner_median_eachnode+Mc;
DmplusMc=Dm+Mc;

subplot(3,3,2)
median_PLR=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,P_LR_median_eachwindow,1,time_flag,timeleg);
xlim([MT_xlim(1),MT_xlim(2)]);
yyaxis left
ylim([MT_ylim(1),MT_ylim(2)]);

subplot(3,3,3)
median_bGR=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,bGR_median_eachnode,2,time_flag,timeleg);
xlim([MT_xlim(1),MT_xlim(2)]);
yyaxis left
ylim([MT_ylim(1),MT_ylim(2)]);
yyaxis right
ylim([b_ylim(1),b_ylim(2)]);

subplot(3,3,5)
median_bTGR=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,bTGR_median_eachnode,3,time_flag,timeleg);
xlim([MT_xlim(1),MT_xlim(2)]);
yyaxis left
ylim([MT_ylim(1),MT_ylim(2)]);
yyaxis right
ylim([b_ylim(1),b_ylim(2)]);

subplot(3,3,6)
median_mcorner=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,mcorner_median_eachnode,4,time_flag,timeleg);
xlim([MT_xlim(1),MT_xlim(2)]);
yyaxis left
ylim([MT_ylim(1),MT_ylim(2)]);
yyaxis right
ylim([m_croner_ylim(1),m_croner_ylim(2)]);

%%
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
%}