function []=CI()
    %% --- Color Settings ---
    color1 = 1/255 * [255, 215, 0; 255, 145, 0; 255, 0, 0];
    color2 = 1/255 * [132, 94, 194; 178, 91, 0; 0, 139, 200];

    %% --- Constants and Parameters ---
    N = 1e3;           % Total number of events
    b = 1;             % b-value of Gutenberg-Richter law
    beta = 2/3 * b;    
    mc_set = 0;        % Magnitude of completeness
    m_corner = 2;      % Corner magnitude
    Mc_set = 10^(1.5 * (mc_set + 6.07)); 
    M_corner = 10^(1.5 * (m_corner + 6.07));

    n0 = 1:N;          % Sample sizes
    dm = 0.01;         % Magnitude step
    m_range = 0:dm:3*log10(N)/b;

    alpha = 90;        % Confidence level (in %)

    %% --- Compute GR Model ---
    M_range_GR = 10.^(1.5 * m_range + 9.1); 
    P_large_GR = 10.^(-b * m_range);

    median_vals_GR = zeros(size(n0));
    lower_bounds_GR = zeros(length(n0), 1);
    upper_bounds_GR = zeros(length(n0), 1);

    for j = 1:length(n0)
        n = n0(j);
        P = binopdf(n, N, P_large_GR);

        [~, idx_max] = max(P);
        median_vals_GR(j) = m_range(idx_max);

        cumulative_P = cumsum(P) / sum(P);
        lower_idx = find(cumulative_P >= (1 - alpha/100)/2, 1, 'first');
        upper_idx = find(cumulative_P >= (1 + alpha/100)/2, 1, 'first');
        lower_bounds_GR(j) = m_range(lower_idx);
        upper_bounds_GR(j) = m_range(upper_idx);
    end

    %% --- Compute TGR Model ---
    M_range_TGR = 10.^(1.5 * (m_range + 6.07));
    P_large_TGR = exp((Mc_set - M_range_TGR)/M_corner) .* (Mc_set ./ M_range_TGR).^beta;

    median_vals_TGR = zeros(size(n0));
    lower_bounds_TGR = zeros(length(n0), 1);
    upper_bounds_TGR = zeros(length(n0), 1);

    for j = 1:length(n0)
        n = n0(j);
        P = binopdf(n, N, P_large_TGR);

        [~, idx_max] = max(P);
        median_vals_TGR(j) = m_range(idx_max);

        cumulative_P = cumsum(P) / sum(P);
        lower_idx = find(cumulative_P >= (1 - alpha/100)/2, 1, 'first');
        upper_idx = find(cumulative_P >= (1 + alpha/100)/2, 1, 'first');
        lower_bounds_TGR(j) = m_range(lower_idx);
        upper_bounds_TGR(j) = m_range(upper_idx);
    end

    %% --- Plot Both Models in One Panel ---
    yy = n0 / N;

    % GR Model
    h(1)=semilogy(median_vals_GR + mc_set, yy, '-', 'Color', [0.2, 0.4, 0.9], 'LineWidth', 2.5);
    hold on;
    fill([lower_bounds_GR; flipud(upper_bounds_GR)] + mc_set, ...
         [yy'; flipud(yy')],[0.2, 0.4, 0.9], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % TGR Model
    [~,ind]=min(abs(median_vals_TGR-m_corner));
    
    h(2)=semilogy(median_vals_TGR + mc_set, yy, '-', 'Color', [0.85, 0.1, 0.1], 'LineWidth', 2.5);
    fill([lower_bounds_TGR; flipud(upper_bounds_TGR)] + mc_set, ...
         [yy'; flipud(yy')], [0.85, 0.1, 0.1], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.3);
    scatter(m_corner,yy(ind),50,[0.85, 0.1, 0.1],'diamond','filled');

    % Formatting
    grid off; box off;
    xlabel('{\itm} - {\itm}_C');
    ylabel('{\itP}(≥{\it m})');    
    set(gca, 'FontSize', 14, 'YScale', 'log');
    ylim([1/N 1]);
    xlim([0 4.5]);
    xticks(0:1:4);
    yticks([1e-3 1e-2 1e-1 1e0]);
  %  legend(h,{'Gutenberg–Richter (GR)','Tapered Gutenberg–Richter (TGR)'}, 'FontSize', 14, 'Location', 'Northeast');
   % text(3,1e-2,'GR','color',[0.2, 0.4, 0.9], 'FontSize', 14);
%     text(1.5, 4e-2, '{\it P}(\geq {\itm}_{max}^{GR or TGR}) = 1/{\itN}', ...
%     'Color', 'k', 'FontSize', 14, 'Interpreter', 'tex');
  %  text(0.4,4e-2,'TGR','color',[0.85, 0.1, 0.1], 'FontSize', 14);
    text(0.35,8e-3,'{\itm}_{corner}','color',[0.85, 0.1, 0.1], 'FontSize', 14);
    legend(h,'GR','TGR','numcolumns',2);
    
end

