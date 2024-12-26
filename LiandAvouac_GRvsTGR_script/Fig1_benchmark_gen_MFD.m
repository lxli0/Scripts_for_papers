clc,clear
color1=1/255*[255,215,0;255,165,0;255,0,0];
color2=1/255*[132 94 194;178 91 0;0 139 200];

b0=[0.6 1 1.5]; % set b-value
L=1e4;
zeta0=0.1:0.1:1;
mmin_set=0;
mmin_true=mmin_set;
Mmin_set=10^(1.5*(mmin_true+6.07));
cc=flip(autumn(length(zeta0)));
fig1=figure('units','normalized','position',[0.1,0.1,0.7,0.5]);
fig2=figure('units','normalized','position',[0.1,0.1,0.7,0.3]);
for i=1:length(b0)
    b=b0(i);
    beta=2/3*b;
    m_range=0:0.02:log10(L)/b+1;
    M_range=10.^(1.5*(m_range+6.07));

    for j=1:length(zeta0)
        zeta=zeta0(j);
        mcorner_TGR=zeta*log10(L)/b;
        M_corner=10^(1.5*(mcorner_TGR+6.07));  
        P_GR=(Mmin_set./M_range).^beta;
        P_TGR=(Mmin_set./M_range).^beta.*exp((Mmin_set-M_range)/M_corner);
        
        figure(fig1);
        subplot(2,length(b0),i)
        semilogy(m_range,L*P_GR,'-k','linewidth',2.5);
        hold on;
        semilogy(m_range,L*P_TGR,'color',cc(j,:),'linewidth',2.5);
        subplot(2,length(b0),i+length(b0))
        semilogy(m_range-mcorner_TGR,L*(P_GR-P_TGR),'color',cc(j,:),'linewidth',2.5);
        hold on;
        
        seq_number=1e2; % number of random catalogs
        Dm_GR_plot=[];
        Dm_TGR_plot=[];
        for k=1:seq_number
            Dm=-1/b.*log10(rand(L,1));
            Dm_GR_plot=[Dm_GR_plot;Dm];
            Dm=gentgr(beta,L,Mmin_set,M_corner)';
            Dm_TGR_plot=[Dm_TGR_plot;Dm];
        end
        figure(fig2);
        subplot(1,length(b0),i)
        Data_GR=plot_GR(Dm_GR_plot,[0 0 0]);hold on;
        Data_TGR=plot_GR(Dm_TGR_plot,cc(j,:));hold on;
    end
    figure(fig1);
    subplot(2,length(b0),i)
    ylim([1 L])
    xlim([0 log10(L)/b])
    set(gca,'fontsize',16)
    grid on;box on;grid minor;
    xlabel('$m$', 'Interpreter', 'latex');
    ylabel('$N(\geq m)$', 'Interpreter', 'latex');
    fortitle = ['\textbf{\textit{b}=', num2str(b), '}'];
    title(fortitle, 'Interpreter', 'latex');
    
    subplot(2,length(b0),i+length(b0))
    ylim([1 L])
   % xlim([0 inf])
    set(gca,'fontsize',16);
    grid on;box on;grid minor;
    xlabel('$m-m_{corner}$', 'Interpreter', 'latex');
    ylabel('$N_{GR}-N_{TGR}$', 'Interpreter', 'latex');
    
    figure(fig2);
    subplot(1,length(b0),i)
    ylim([1e-4 1])
    xlim([0 log10(L)/b])
    set(gca,'fontsize',16)
    grid on;box on;grid minor;
    xlabel('$m$', 'Interpreter', 'latex');
    ylabel('$S = P(\geq m)$', 'Interpreter', 'latex');
    fortitle = ['\textbf{\textit{b}=', num2str(b), '}'];
    title(fortitle, 'Interpreter', 'latex');
    
end
figure(fig1);
forleg = arrayfun(@(x) ['$\zeta=' num2str(x) '$'], zeta0, 'UniformOutput', false);
legend(forleg, 'NumColumns', 2, 'Location', 'Northeast', 'Interpreter', 'latex');

function [Dm]=gentgr(beta,L,Mmin,M_corner)
        U10=rand(1,L);
        U20=rand(1,L);
        U1=Mmin-M_corner*log(U10);
        U2=Mmin*U20.^(-1/beta);
        DM=min(U1,U2);
        Dm=2/3*log10(DM)-6.07;    
end

function [pl]=plot_GR(Dm,color)
    L=length(Dm);
    dm=0.1;
    m=-0.5:dm:max(Dm);
    n0=hist(Dm,m);  
    cn0=[];
    cn0(1)=L;
    for i=2:length(m)
        cn0(i)=L-sum(n0(1:i-1)); 
    end
    pl=semilogy(m,cn0/L,'color',color,'linewidth',2.5);
    hold on;
end