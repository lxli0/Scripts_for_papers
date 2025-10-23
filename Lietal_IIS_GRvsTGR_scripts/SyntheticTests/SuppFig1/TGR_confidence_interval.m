clc,clear
close all
color1=1/255*[255,215,0;255,165,0;255,0,0;];
color2=1/255*[132 94 194;178 91 0;0 139 200];
N=1e3;
b=1;
beta=2/3*b;
mc_set=0;
m_corner=3;
Mc_set=10^(1.5*(mc_set+6.07));
M_corner=10^(1.5*(m_corner+6.07));    
n0=1:1:N;
dm=0.01;
m_range=0:dm:3*log10(N)/b;
M_range=10.^(1.5*(m_range+6.07)); 
P_large=exp((Mc_set-M_range)/M_corner).*(Mc_set./M_range).^beta;
alpha=[90,95,99];
for j=1:length(n0)
    n=n0(j);
    P=binopdf(n,N,P_large);
    [z0,h0]=max(P);
    meadian(j)=m_range(h0);
    TT=cumsum(P)/sum(P);
    for i=1:length(alpha)
        [z1,h1]=min((abs(TT-(1-alpha(i)/100)/2)));
        [z2,h2]=min((abs(TT-(1+alpha(i)/100)/2)));
        small(j,i)=m_range(h1);
        large(j,i)=m_range(h2);
    end
end
%%
figure('units','normalized','position',[0.1,0.1,0.25,0.25])
L=1e3;% number of random catalog
for i=1:L
    Dm=gentgr(beta,L,Mc_set,M_corner)-0.05;
    plotgr(Dm,[0.75 0.75 0.75],1);
end
hold on;
for i=1:length(alpha)
    GR_con(i)=semilogy(small(:,i)+mc_set,n0,'-','color',color1(i,:),'linewidth',2.5);
    hold on;
    semilogy(large(:,i)+mc_set,n0,'-','color',color1(i,:),'linewidth',2.5);
    hold on;
end
ME=semilogy(meadian+mc_set,n0,'-','color',color2(3,:),'linewidth',2.5);
grid on;box on;
xlabel('Magnitude');
ylabel('Cumulative Number');
set(gca,'fontsize',16);
ylim([1 inf]);
xlim([0 2*log10(N)/b])
%legend([ME,GR_con],'Mode','90%','95%','99%')
hold on;
xticks([0:1:6]);
title('TGR, {\itb} = 1, {\itm}_{corner} = 3');



%%
function [pl]=plotgr(Dm,color,linewidth)
    Catalogsize=length(Dm);
    dm=0.1;
    m=0:dm:max(Dm);
    n0=hist(Dm,m); 
    cn0=[];
    cn0(1)=Catalogsize;
    for i=2:length(m)
        cn0(i)=Catalogsize-sum(n0(1:i-1)); 
    end
    pl=semilogy(m,cn0,'color',color,'linewidth',linewidth);
    hold on;
end


function [Dm]=gentgr(beta,L,Mmin,M_corner)
        U10=rand(1,L);
        U20=rand(1,L);
        U1=Mmin-M_corner*log(U10);
        U2=Mmin*U20.^(-1/beta);
        DM=min(U1,U2);
        Dm=2/3*log10(DM)-6.07;    
end
