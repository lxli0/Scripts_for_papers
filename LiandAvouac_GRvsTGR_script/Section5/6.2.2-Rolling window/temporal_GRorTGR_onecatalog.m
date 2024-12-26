figure('units','normalized','position',[0.1,0.1,0.5,0.4])
clc,clear
Nmax=3000;
mmin=0;

x =-1/0.9.*log10(rand(1,Nmax/3));

beta=2/3;
Mmin=10^(1.5*(mmin+6.07));
mcorner_TGR=2;
M_corner=10^(1.5*(mcorner_TGR+6.07));

U10=rand(1,Nmax/3);
U20=rand(1,Nmax/3);
U1=Mmin-M_corner*log(U10);
U2=Mmin*U20.^(-1/beta);
DM=min(U1,U2);
x=[x,2/3*log10(DM)-6.07];                   

x =[x, -1/1.1*log10(rand(1,Nmax/3))];
yyaxis left
stem(1:Nmax,x,'k');
ylabel('Magntidue');
ax = gca; 
ax.YColor = 'k';
ylim([0 5])

yyaxis right
color1=1/255*[255,215,0;255,165,0;255,0,0;];
color2=1/255*[132 94 194;178 91 0;0 139 200];
x_values = 1:1:Nmax/3;
y_values = ones(1, length(x_values));
fill([x_values fliplr(x_values)], [zeros(1,length(x_values)) fliplr(y_values)], color1(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on;
text(Nmax/6, 0.95, 'GR', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
text(Nmax/6, 0.88, '$b$=0.9', 'HorizontalAlignment', 'center', 'FontSize', 16, 'Interpreter', 'latex', 'FontWeight', 'bold');
x_values = Nmax/3+1:1:Nmax/3*2;
y_values = ones(1, length(x_values));
fill([x_values fliplr(x_values)], [zeros(1,length(x_values)) fliplr(y_values)], color2(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on;
text(Nmax/2, 0.95, 'TGR', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
text(Nmax/2, 0.88, '$b$=1.0', 'HorizontalAlignment', 'center', 'FontSize', 16, 'Interpreter', 'latex', 'FontWeight', 'bold');
text(Nmax/2, 0.83, '$m_{corner}$=2.5', 'HorizontalAlignment', 'center', 'FontSize', 16, 'Interpreter', 'latex', 'FontWeight', 'bold');
x_values = Nmax/3*2+1:1:Nmax;
y_values = ones(1, length(x_values));
fill([x_values fliplr(x_values)], [zeros(1,length(x_values)) fliplr(y_values)], color1(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on;
text(Nmax/6*5, 0.95, 'GR', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
text(Nmax/6*5, 0.88, '$b$=1.1', 'HorizontalAlignment', 'center', 'FontSize', 16, 'Interpreter', 'latex', 'FontWeight', 'bold');

Ntotal=length(x);
min_size=100;
max_size=1e3;
w0=ceil(Ntotal/max_size)+1:1:min([ceil(Ntotal/min_size)-1,ceil(Ntotal/max_size)+1+4]); % number of windows

[P_LR_median_eachwindow]=get_temporal_MFD(x,mmin,min_size,max_size,w0);
col=autumn(length(w0));
for k=1:length(w0)
        h(k)=plot(1:Nmax,P_LR_median_eachwindow(k,:),'-','color',col(k,:),'linewidth',2.5);hold on;
        leg{k}=['$w=$',num2str(w0(k))];
end
meadian_plt_variable=median(P_LR_median_eachwindow,1);
h(k+1)=plot(1:Nmax,meadian_plt_variable,'-b','linewidth',2.5);hold on;
leg{k+1}=['median'];
xlabel('# of events')
ylabel('$p$','Interpreter','latex')
set(gca,'fontsize',16);
grid on;box on;grid minor;
ylim([0 1]);
ax = gca; 
ax.YColor = 'k';
legend(h,leg, 'Interpreter', 'latex','location','northwest','Numcolumns',3);

%%

function [P_LR_median_eachwindow]=get_temporal_MFD(x,mmin,min_size,max_size,w0)
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