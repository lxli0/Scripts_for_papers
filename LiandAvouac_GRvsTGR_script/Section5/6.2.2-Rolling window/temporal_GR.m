figure('units','normalized','position',[0.1,0.1,0.6,0.3])
%{
clc,clear
Nmax=2000;
mmin=0;
for sk=1:100
    sk
    x =-1/0.9.*log10(rand(1,Nmax/2));
    x =[x, -1/1.1*log10(rand(1,Nmax/2))];
    Ntotal=length(x);
    min_size=100;
    max_size=1e3;
    w0=ceil(Ntotal/max_size)+1:1:min([ceil(Ntotal/min_size)-1,ceil(Ntotal/max_size)+1+4]); % number of windows
    b_median_node(sk,:)=get_temporal_GGR(x,mmin,min_size,max_size,w0);
end
save('temp_GR_1.mat')
%}
clc,clear
load('temp_GR_1.mat');
subplot(1,2,1)
color1=1/255*[255,215,0;255,165,0;255,0,0;];
color2=1/255*[132 94 194;178 91 0;0 139 200];

for sk=1:100
    h=plot(b_median_node(sk,:),'color',[0.75 0.75 0.75],'linewidth',2);hold on;
end
hh=plot(mean(b_median_node,1),'k','linewidth',2.5);hold on;
hhh=plot([1 1e3 1e3 2e3],[0.9 0.9 1.1 1.1],'color',color1(2,:),'linewidth',2);
xlabel('# of events')
ylabel('$b$','Interpreter','latex')
set(gca,'fontsize',16);
grid on;box on;grid minor;
ylim([0.8 1.2]);
legend([h,hh,hhh],'Random catalog','Mean','True','location','southeast');

%%
clc,clear
Nmax=2000;
mmin=0;
for sk=1:100
    b_true=0.9+0.2*(1:Nmax)/Nmax;
    x =-1./b_true.*log10(rand(1,Nmax));
    Ntotal=length(x);
    min_size=100;
    max_size=1e3;
    w0=ceil(Ntotal/max_size)+1:1:min([ceil(Ntotal/min_size)-1,ceil(Ntotal/max_size)+1+4]); % number of windows
    b_median_node(sk,:)=get_temporal_GGR(x,mmin,min_size,max_size,w0);
end
save('temp_GR_2.mat')
%}
clc,clear
load('temp_GR_2.mat');
subplot(1,2,2)
color1=1/255*[255,215,0;255,165,0;255,0,0;];
color2=1/255*[132 94 194;178 91 0;0 139 200];

for sk=1:100
    h=plot(b_median_node(sk,:),'color',[0.75 0.75 0.75],'linewidth',2);hold on;
end
hh=plot(mean(b_median_node,1),'k','linewidth',2.5);hold on;
b_true=0.9+0.2*(1:Nmax)/Nmax;
hhh=plot(1:Nmax,b_true,'color',color1(2,:),'linewidth',2);
xlabel('# of events')
ylabel('$b$','Interpreter','latex')
set(gca,'fontsize',16);
grid on;box on;grid minor;
ylim([0.8 1.2]);

%%

function [bGR_median]=get_temporal_GGR(x,mmin,min_size,max_size,w0)
    Ntotal=length(x);
    num_rand=100; % for each w0, number of random partition
    bGR_median_eachnode=zeros(length(w0),Ntotal);
    
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
                bGR(i,node(i,j):node(i,j+1))=b_GR; 
            end
        end
        bGR_median_eachnode(k,:)=median(bGR,1);      
    end
    bGR_median=median(bGR_median_eachnode,1);
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
