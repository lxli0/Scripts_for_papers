clc,clear
color1=1/255*[255,215,0;255,165,0;255,0,0;];
color2=1/255*[132 94 194;178 91 0;0 139 200];
N=1e3;
b=1;
beta=2/3*b;
mc_set=0;
m_corner=2;
Mc_set=10^(1.5*(mc_set+6.07));
M_corner=10^(1.5*(m_corner+6.07));    
n0=1:1:N;
dm=0.01;
m_range=0:dm:3*log10(N)/b;
M_range=10.^(1.5*m_range+9.1); 
P_large=10.^(-b*m_range);
alpha=[90,95,99];
L=1e5;
%%
figure('units','normalized','position',[0.1,0.1,0.3,0.3])
Mag0=zeros(N,L);
for i=1:L
    Dm=gentgr(beta,L,Mc_set,M_corner)-0.05;% generate random magnitude
    dDm=diff(Dm);
    dDm_pos=dDm(dDm>=-0.05);
    dDm_pos=dDm_pos(1:N);
    if i<=1e3
        plotgr(dDm_pos,[0.75 0.75 0.75],1);
    end
    
    Dm_cal=roundn(dDm_pos+0.05,log10(0.1));
    Dm_cal=sort(Dm_cal, 'descend');
    Mag0(:,i)=Dm_cal;
    max_M(i)=max(Dm_cal);
end

for j=1:N
    Mag=sort(Mag0(j,:),'descend');    
    meadian(j)=mode(Mag);
    for i=1:length(alpha)
        small(j,i)=Mag(floor(alpha(i)/100*L));
        large(j,i)=Mag(floor((1-alpha(i)/100)*L));
    end
end

for i=1:length(alpha)
    GR_con(i)=semilogy(small(:,i)+mc_set,n0,'-','color',color1(i,:),'linewidth',2.5);
    hold on;
    semilogy(large(:,i)+mc_set,n0,'-','color',color1(i,:),'linewidth',2.5);
    hold on;
end
ME=semilogy(meadian+mc_set,n0,'-','color',color2(3,:),'linewidth',2.5);
grid on;box on;
xlabel('Magnitude');
ylabel('Culmulative Number');
set(gca,'fontsize',16);
ylim([1 inf]);
xlim([0 2*log10(N)/b])
%legend([ME,GR_con],'Mode','90%','95%','99%')
hold on;

grid on;box on;
xlabel('Magnitude Difference');
ylabel('Culmulative Number');
set(gca,'fontsize',16);
ylim([1 inf]);
xlim([0 2*log10(N)/b])

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
