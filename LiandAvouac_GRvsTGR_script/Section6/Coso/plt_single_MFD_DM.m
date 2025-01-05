% M is initial catalog, Dm is the adjuasted catalog (above mc, and minus
% mc)
function []=plt_single_MFD_DM(M,Dm,mc,b_GR)
    %% Plot cum- and non-culm distribution
    L=length(M);
    dm=0.1;
    m=floor(min(M)/dm)*dm:dm:max(M);
    n0=hist(M,m);   
    cn0(1)=L;
    for i=2:length(m)
        cn0(i)=length(M(M>=m(i)-1e-12)); 
    end
    n=log10(n0);
    cn=log10(cn0);
    color1=[0.75 0.75 0.75; 0.25 0.25 0.25];
    color2=1/255*[132 94 194;178 91 0;0 139 200];
    plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
                'markersize', 5, 'markerfacecolor', color1(1,:), ...
                'MarkerEdgeColor', color1(1,:));
    hold on;
    plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
                'markersize', 5, 'markerfacecolor', color1(2,:), ...
                'MarkerEdgeColor', color1(2,:));
    hold on;
    ylimmax=ceil(L/10^floor(log10(L)))*10^floor(log10(L));
    ylim([1 2e4]); 
    yticks([1 1e1 1e2 1e3 1e4]);
    if isnan(mc)==0
        %% Plot the best fit
        m=0:dm:max(Dm);
        Moment_min=10^(1.5*(min(m)+6.07));
        F_GR=(Moment_min./10.^(3/2*(m+6.07))).^(2/3*b_GR);
        fit_GR=semilogy(m+mc,length(Dm)*F_GR,'color',color2(2,:),'linewidth',2.5);
        hold on;

        scatter(mc,length(Dm)*F_GR(1)*1.5,100,1/255*[255,165,0],'v','filled');
        hold on;
        %% Plot the confidential interval
        GR_con=GR_confidential_interval(b_GR,length(Dm),mc);
       %{
        if plt_num==1
            legend([plot_n0,plot_cn0,fit_GR,fit_TGR,GR_con],'Non-cumulative','Cumulative','Fit_G_R', ...
                'Fit_T_G_R','90% confidential interval','Location','Northeast','NumColumns',1);
        end
%}
        %%
        xlim([0 5]);
        xticks([0 1 2 3 4 5]);
        ax = gca; 
        xTextPos = 0.95*ax.XLim(2); 
        yTextPos = 0.9*10^(log10(ax.YLim(2)));
        
        textStr = sprintf('$m_c''$: %.1f\\\\\n$b_{GR}$: %.2f', mc, b_GR);
        text(xTextPos, yTextPos, textStr,'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 16,'Interpreter', 'latex');
    end
    grid on;box on;grid minor;
    xlabel('Magnitude Difference');
    ylabel('Number');
    set(gca,'fontsize',16);
    hold on;   
end

function [pl]=GR_confidential_interval(b,N,mc_set)
    color1=1/255*[132 94 194;178 91 0;0 139 200];
    n0=1:1:N;
    dm=0.01;
    m_range=0:dm:3*log10(N)/b;
    P_large=10.^(-b*m_range);
    alpha=90;
    for j=1:length(n0)
        n=n0(j);
        P=binopdf(n,N,P_large);
        [z0,h0]=max(P);
        meadian(j)=m_range(h0);
        TT=cumsum(P)/sum(P);
        [z1,h1]=min((abs(TT-(1-alpha/100)/2)));
        [z2,h2]=min((abs(TT-(1+alpha/100)/2)));
        small(j)=m_range(h1);
        large(j)=m_range(h2);
    end
    pl=semilogy(small+mc_set,n0,'--','color',color1(2,:),'linewidth',1.5);
    hold on;
    semilogy(large+mc_set,n0,'--','color',color1(2,:),'linewidth',1.5);
    hold on;
end


