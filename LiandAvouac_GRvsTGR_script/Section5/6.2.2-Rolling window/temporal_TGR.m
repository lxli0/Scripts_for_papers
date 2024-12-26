figure('units','normalized','position',[0.1,0.1,0.6,0.6])
%{
clc,clear
Nmax=2000;
mmin=0;
Mmin=10^(1.5*(mmin+6.07));
beta=2/3;
parpool(8);
parfor sk=1:100
    sk
    
    mcorner_TGR=1.5;
    M_corner=10^(1.5*(mcorner_TGR+6.07));    
    [Dm]=gentgr(beta,Nmax/2,Mmin,M_corner)
    x=Dm;
    
    mcorner_TGR=2;
    M_corner=10^(1.5*(mcorner_TGR+6.07));    
    [Dm]=gentgr(beta,Nmax/2,Mmin,M_corner)
    x=[x,Dm];
   
    Ntotal=length(x);
    min_size=100;
    max_size=1e3;
    w0=ceil(Ntotal/max_size)+1:1:min([ceil(Ntotal/min_size)-1,ceil(Ntotal/max_size)+1+4]); % number of windows

    [P_LR_median(sk,:),bGR_median(sk,:),bTGR_median(sk,:),mcorner_median(sk,:)]=get_temporal_MFD(x,mmin,min_size,max_size,w0)

end
delete(gcp('nocreate'))
save('temp_TGR_1.mat')
%}
clc,clear
load('temp_TGR_1.mat');
color1=1/255*[255,215,0;255,165,0;255,0,0;];
color2=1/255*[132 94 194;178 91 0;0 139 200];
subplot(2,2,1)
for sk=1:100
    h=plot(bTGR_median(sk,:),'color',[0.75 0.75 0.75],'linewidth',2);hold on;
end
hh=plot(mean(bTGR_median,1),'k','linewidth',2.5);hold on;
hhh=plot([1 1e3 1e3 2e3],[1 1 1 1],'color',color1(2,:),'linewidth',2);
xlabel('# of events')
ylabel('$b$','Interpreter','latex')
set(gca,'fontsize',16);
grid on;box on;grid minor;
ylim([0.8 1.2]);
legend([h,hh,hhh],'Random catalog','Mean','True','location','southeast');

subplot(2,2,3)
for sk=1:100
    h=plot(mcorner_median(sk,:),'color',[0.75 0.75 0.75],'linewidth',2);hold on;
end
hh=plot(mean(mcorner_median,1),'k','linewidth',2.5);hold on;
hhh=plot([1 1e3 1e3 2e3],[1.5 1.5 2 2],'color',color1(2,:),'linewidth',2);
xlabel('# of events')
ylabel('$m_{corner}$','Interpreter','latex')
set(gca,'fontsize',16);
grid on;box on;grid minor;
ylim([1 2.5]);

%%
%{
clc,clear
Nmax=2000;
mmin=0;
Mmin=10^(1.5*(mmin+6.07));
beta=2/3;
parpool(8);
parfor sk=1:100
    sk
    mcorner_TGR=zeros(Nmax,1);
    for i = 1:Nmax
        mcorner_TGR(i)=1.5+0.5*i/Nmax;
    end

    x = zeros(Nmax, 1);
    for i = 1:Nmax
        M_corner=10^(1.5*(mcorner_TGR(i)+6.07));                   
        x(i)=gentgr(beta,1,Mmin,M_corner);
    end
   
    Ntotal=length(x);
    min_size=100;
    max_size=1e3;
    w0=ceil(Ntotal/max_size)+1:1:min([ceil(Ntotal/min_size)-1,ceil(Ntotal/max_size)+1+4]); % number of windows

    [P_LR_median(sk,:),bGR_median(sk,:),bTGR_median(sk,:),mcorner_median(sk,:)]=get_temporal_MFD(x,mmin,min_size,max_size,w0)

end
delete(gcp('nocreate'))
save('temp_TGR_2.mat')
%}

clc,clear
load('temp_TGR_2.mat');
color1=1/255*[255,215,0;255,165,0;255,0,0;];
color2=1/255*[132 94 194;178 91 0;0 139 200];
subplot(2,2,2)
for sk=1:100
    h=plot(bTGR_median(sk,:),'color',[0.75 0.75 0.75],'linewidth',2);hold on;
end
hh=plot(mean(bTGR_median,1),'k','linewidth',2.5);hold on;
hhh=plot([1 1e3 1e3 2e3],[1 1 1 1],'color',color1(2,:),'linewidth',2);
xlabel('# of events')
ylabel('$b$','Interpreter','latex')
set(gca,'fontsize',16);
grid on;box on;grid minor;
ylim([0.8 1.2]);

subplot(2,2,4)
for sk=1:100
    h=plot(mcorner_median(sk,:),'color',[0.75 0.75 0.75],'linewidth',2);hold on;
end
hh=plot(mean(mcorner_median,1),'k','linewidth',2.5);hold on;
mcorner_TGR=zeros(Nmax,1);
for i = 1:Nmax
    mcorner_TGR(i)=1.5+0.5*i/Nmax;
end
hhh=plot(1:Nmax,mcorner_TGR,'color',color1(2,:),'linewidth',2);
xlabel('# of events')
ylabel('$m_{corner}$','Interpreter','latex')
set(gca,'fontsize',16);
grid on;box on;grid minor;
ylim([1 2.5]);
%}
%%

function [P_LR_median,bGR_median,bTGR_median,mcorner_median]=get_temporal_MFD(x,mmin,min_size,max_size,w0)
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
                [b_GR,loglikelihood_GR]=Estimation_GR(xx,mmin);
                [b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR_gridsearch_continuous(xx,mmin);
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
    P_LR_median=median(P_LR_median_eachwindow,1);
    bGR_median=median(bGR_median_eachnode,1);
    bTGR_median=median(bTGR_median_eachnode,1);
    mcorner_median=median(mcorner_median_eachnode,1);
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

function [Dm]=gentgr(beta,L,Mmin,M_corner)
        U10=rand(1,L);
        U20=rand(1,L);
        U1=Mmin-M_corner*log(U10);
        U2=Mmin*U20.^(-1/beta);
        DM=min(U1,U2);
        Dm=2/3*log10(DM)-6.07;    
end