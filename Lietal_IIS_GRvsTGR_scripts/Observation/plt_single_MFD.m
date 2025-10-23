% M is initial catalog, Dm is the adjuasted catalog (above mc, and minus
% mc)
function []=plt_single_MFD_for_fig1(M,Dm,mc,b_GR,b_TGR,Mcorner_TGR_minus_mc,P_LR,alpha_GR,alpha_TGR,leg_flag)
    %% Plot cum- and non-culm distribution
    L=length(M);
    dm=0.1;
    m=floor(min(M)/dm)*dm:dm:max(M);
    if max(m)<max(M)
        m=[m,max(M)];
    end
    n0=hist(M,m);   
    for i=1:length(m)
        cn0(i)=length(find((M-m(i))>=-1e-10));
    end
    n=log10(n0);
    cn=log10(cn0);
    color1=[0.75 0.75 0.75; 0.25 0.25 0.25];
    color2=1/255*[132 94 194;178 91 0;0 139 200];
    plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
                'markersize', 6.5, 'markerfacecolor', color1(1,:), ...
                'MarkerEdgeColor', color1(1,:));
    hold on;
    plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
                'markersize', 6.5, 'markerfacecolor', color1(2,:), ...
                'MarkerEdgeColor', color1(2,:));
    hold on;
    
    if isnan(mc)==0
        %% Plot the best fit
        m=0:dm:max(Dm);
        if max(m)<max(Dm)
            m=[m,max(Dm)];
        end
        Moment_min=10^(1.5*(min(m)+6.07));
        Moment_corner=10^(1.5*(Mcorner_TGR_minus_mc+6.07));
        F_GR=(Moment_min./10.^(3/2*(m+6.07))).^(2/3*b_GR);
        F_TGR=(Moment_min./10.^(3/2*(m+6.07))).^(2/3*b_TGR).*exp((Moment_min-10.^(3/2*(m+6.07)))/ Moment_corner);
        
         %% Plot the confidential interval
        GR_con=GR_confidential_interval(b_GR,length(Dm),mc);

        
       % fit_GR=semilogy(m+mc,length(Dm)*F_GR,'color',color2(3,:),'linewidth',3);
        hold on;
        fit_TGR=semilogy(m+mc,length(Dm)*F_TGR,'--','color',[0.85, 0.1, 0.1],'linewidth',3);
        hold on;
        scatter(mc,length(Dm)*F_GR(1)*1.5,100,1/255*[255,165,0],'v','filled');
        hold on;
    
       %     legend([plot_n0,plot_cn0,fit_GR,fit_TGR,GR_con],'Non-cumulative','Cumulative','Fit (GR)', ...
       %         'Fit (T_G_R)','90%','Location','Northeast','NumColumns',1);
 
        %%
        % Set x-axis limits and ticks
        xMin = round(min(M));
        xMax = ceil(max(M));
        xlim([xMin xMax]);
        xticks(xMin:xMax); % Sets ticks at every integer between xMin and xMax
        xticklabels(arrayfun(@num2str, xMin:xMax, 'UniformOutput', false));

        
        % Set y-axis limits and ticks
        yMin = 0; % Assuming the minimum y value starts from 1
        yMax = ceil(log10(L));  % L is the maximum value in the original y data
        ylim([10^yMin, 10^yMax]);

        % Set y-axis ticks and labels
        yticks(10.^(yMin:yMax)); % Logarithmic ticks
        yticklabels(arrayfun(@(y) sprintf('10^{%d}', round(log10(y))), yticks, 'UniformOutput', false));

        
        ax = gca; 
        if ax.XLim(2)>0
            xTextPos = 0.95*ax.XLim(2); 
        else
            xTextPos = 1.05*ax.XLim(2); 
        end
        yTextPos = 0.9*10^(log10(ax.YLim(2)));
if leg_flag==1
textStr = sprintf('$m_{\\mathrm{C}}$: %.1f, $b_{\\mathrm{GR}}$: %.2f, $b_{\\mathrm{TGR}}$: %.2f\\\\\n$m_{\\mathrm{corner}}$: %.1f, $p_{\\mathrm{LRT}}$: %.2f\\\\\n$P(p_{\\mathrm{LRT}}^{\\mathrm{GR}} < p_{\\mathrm{LRT}}^{\\mathrm{OBS}})$: %.2f\\\\\n$P(p_{\\mathrm{LRT}}^{\\mathrm{TGR}} > p_{\\mathrm{LRT}}^{\\mathrm{OBS}})$: %.2f', ...
    mc, b_GR, b_TGR, Mcorner_TGR_minus_mc+mc, P_LR, alpha_GR, alpha_TGR);
      
        text(xTextPos, yTextPos, textStr,'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 13,'Interpreter', 'latex');
end

    end
    grid on;box off;
        xlabel('\it m');
        ylabel('\it N');
    set(gca,'fontsize',13);
    hold on;   
end
function [pl] = GR_confidential_interval(b, N, mc_set)
    color1 = 1/255 * [132 94 194; 178 91 0; 0 139 200];
    n0 = 1:1:N;
    dm = 0.01;
    m_range = 0:dm:3*log10(N)/b;
    P_large = 10.^(-b * m_range);
    alpha = 90;

    for j = 1:length(n0)
        n = n0(j);
        P = binopdf(n, N, P_large);
        [~, h0] = max(P);
        meadian(j) = m_range(h0);
        TT = cumsum(P) / sum(P);
        [~, h1] = min(abs(TT - (1 - alpha/100)/2));
        [~, h2] = min(abs(TT - (1 + alpha/100)/2));
        small(j) = m_range(h1);
        large(j) = m_range(h2);
    end

    % Shade the confidence interval
    x_fill = [small + mc_set, fliplr(large + mc_set)];
    y_fill = [n0, fliplr(n0)];
    fill(x_fill, y_fill, color1(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    hold on;

    % Optionally plot the median line
    pl = plot(meadian + mc_set, n0, '-', 'Color',[0.2, 0.4, 0.9], 'LineWidth', 3);
    hold on;
end



