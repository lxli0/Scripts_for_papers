
clc;
clear;
close all;

% === Load file names ===
namelist = dir('./Data/*.txt');
file_name = {namelist.name};
case_num = length(file_name);
sites = cell(case_num, 1);

for i = 1:case_num
    [~, sites{i}] = fileparts(file_name{i});
end

% === Load summary and custom colormap ===
load('Case_Summary.mat');
custom_cmap = get_colormap();
Data_Flag = load('./Data/_Data_Flag.md');
Depth_data = load('./Data/_Depth.md');

[~, sort_idx] = sort(Case_Summary(:,12), 'descend');
sites = sites(sort_idx);
Case_Summary = Case_Summary(sort_idx,:);
Data_Flag = Data_Flag(sort_idx,:);
Depth_data = Depth_data(sort_idx,:);

% === Create figures ===
fig1 = figure('units', 'normalized', 'position', [0.05, 0.05, 0.4, 0.8]);

Numrow = 6;
Numcol = 4;
quantity = zeros(case_num, Numrow*Numcol); % Preallocate quantity matrix

% === Loop over cases ===
for i = 1:case_num
    site = sites{i};
    file = ['./Data/', site, '.txt'];
    D = load(file); 
    D = sortrows(D, 1);

    mc = Case_Summary(i,1);
    Num = Case_Summary(i,2);
    b_GR = Case_Summary(i,3);
    b_TGR = Case_Summary(i,6);
    m_corner = Case_Summary(i,9);
    P_LR = Case_Summary(i,12);
    Obs_max = Case_Summary(i,15);
    mode_max_GR = Case_Summary(i,16);
    mode_max_TGR = Case_Summary(i,19);
    a_GR = log10(Num) + b_GR * mc;

    D = D(D(:,5) >= mc-1e-12, :);
    T0 = D(:,1);
    X0 = D(:,2);
    X0 = X0 - mean(X0);
    Y0 = D(:,3);
    Y0 = Y0 - mean(Y0);
    Z0 = D(:,4);
    Z0 = Z0 - mean(Z0);
    M0 = D(:,5);
    
    Inje_V = Data_Flag(i,1);
    peak_Q = Data_Flag(i,2);
    peak_P = Data_Flag(i,3);
    hori_flag = Data_Flag(i,4);
    depth_flag = Data_Flag(i,5);
    time_flag = Data_Flag(i,6);
    mag_flag = Data_Flag(i,7);
    sci_flag = Data_Flag(i,9);
    
    if hori_flag == 1
        X0 = X0 * 110e3;
        Y0 = Y0 * 110e3;
    end
    if depth_flag == 1
        Z0 = Z0 * 1e3;
    end
    Z0 = -abs(Z0);
    if time_flag == 2
        T0 = (T0 - min(T0)) * 365.25;
    end
    
    duration = max(T0) - min(T0);
    diff_T = diff(T0);
    nonzero_vals = diff_T(diff_T ~= 0);
    min_nonzero = min(abs(nonzero_vals));
    diff_T(diff_T == 0) = min_nonzero / 2;
    diff_M = diff(M0);
    
    CV_dT = std(diff_T) / mean(diff_T);
    num_diff = length(diff_T);
    LV_dT = 0;
    x_diffT = 0;
    x2_diffT = 0;
    x_diffM = 0;
    x2_diffM = 0;
    
    for jk = 1:num_diff-1
        LV_dT = LV_dT + 3*(diff_T(jk) - diff_T(jk+1))^2 / (diff_T(jk) + diff_T(jk+1))^2;
        x_diffT = x_diffT + diff_T(jk+1) - diff_T(jk);
        x2_diffT = x2_diffT + (diff_T(jk+1) - diff_T(jk))^2;
        x_diffM = x_diffM + diff_M(jk+1) - diff_M(jk);
        x2_diffM = x2_diffM + (diff_M(jk+1) - diff_M(jk))^2;
    end
    
    LV_dT = LV_dT / (num_diff-1);
    x_diffT = x_diffT / (num_diff-1);
    x2_diffT = x2_diffT / (num_diff-1);
    x_diffM = x_diffM / (num_diff-1);
    x2_diffM = x2_diffM / (num_diff-1);
    
    if ~isnan(X0)
        [Lmax, Lmedian, Lmin, Vellipsoid, coverage] = get_spatial_ellipsoid(X0, Y0, Z0, 2);
        df = fractdime(X0, Y0, Z0, 0);
    else
        Lmax = nan;
        Lmedian = nan;
        Lmin = nan;
        Vellipsoid = nan;
        df = nan;
    end
    
    [~, ind] = max(M0);
    location_mmax = ind / length(M0);
    
    % Plotting in fig1
    color = custom_cmap(1 + floor(round(P_LR, 2) * 100), :);
    figure(fig1);
    
    % Subplot 1: Injection Volume
    t = 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, log10(Inje_V), color, num2str(i), sci_flag);
    quantity(i, t) = log10(Inje_V);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('log_{10}{\itV}_{inj} (m^3)');
    
    % Subplot 2: Peak Flow Rate
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, log10(peak_Q), color, num2str(i), sci_flag);
    quantity(i, t) = log10(peak_Q);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('log_{10}{\itQ}_{max} (l/s)');
    
    % Subplot 3: Peak Pressure
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, peak_P, color, num2str(i), sci_flag);
    quantity(i, t) = peak_P;
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('{\itP}_{max} (MPa)');
    
    % Subplot 4: Event Count
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, log10(Num), color, num2str(i), sci_flag);
    quantity(i, t) = log10(Num);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('log_{10}[{\itN}(≥{\itm}_C)]');
    
    % Subplot 5: GR a-value
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, a_GR, color, num2str(i), sci_flag);
    quantity(i, t) = a_GR;
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('{\ita}_{GR}');
    
    % Subplot 6: GR b-value
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, b_GR, color, num2str(i), sci_flag);
    quantity(i, t) = b_GR;
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('{\itb}_{GR}');
    
    % Subplot 7: TGR b-value
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, b_TGR, color, num2str(i), sci_flag);
    quantity(i, t) = b_TGR;
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('{\itb}_{TGR}');
    
    % Subplot 8: Corner Magnitude
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, m_corner, color, num2str(i), sci_flag);
    quantity(i, t) = m_corner;
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('{\itm}_{corner}');
    
%     % Subplot 9: a/b - m_corner
%     t = t + 1;
%     subplot(Numrow, Numcol, t); hold on;
%     plt_scatter(P_LR, a_GR/b_GR - m_corner, color, num2str(i), sci_flag);
%     quantity(i, t) = a_GR/b_GR - m_corner;
%     set(gca, 'fontsize', 12);
%     xlabel('{\itp}_{LRT}'); ylabel('{\ita}_{GR}/{\itb}_{GR}-{\itm}_{corner}');
%     
    % Subplot 10: Maximum Observed Magnitude
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, Obs_max, color, num2str(i), sci_flag);
    quantity(i, t) = Obs_max;
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('{\itm}_{max}');
    
    % Subplot 11: Duration
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, log10(duration), color, num2str(i), sci_flag);
    quantity(i, t) = log10(duration);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('log_{10}Duration (d)');
    
    % Subplot 12: Normalized a-value
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    normalized_aGR = log10(10^a_GR/duration);
    plt_scatter(P_LR, normalized_aGR, color, num2str(i), sci_flag);
    quantity(i, t) = normalized_aGR;
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('{\ita} (d^{-1})');
    
    % Subplot 13: Seismogenic Index
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    SI = a_GR - log10(Inje_V);
    plt_scatter(P_LR, SI, color, num2str(i), sci_flag);
    quantity(i, t) = SI;
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('\Sigma');
    
    % Subplot 14: CV of Time Differences
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, log10(CV_dT), color, num2str(i), sci_flag);
    quantity(i, t) = log10(CV_dT);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('log_{10}CV_{\Delta\itt}');
    
    % Subplot 15: LV of Time Differences
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, log10(LV_dT), color, num2str(i), sci_flag);
    quantity(i, t) = log10(LV_dT);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('log_{10}LV_{\Delta\itt}');
    
    % Subplot 16: Median Depth
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, Depth_data(i,1), color, num2str(i), sci_flag);
    quantity(i, t) = Depth_data(i,1);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('Depth_{median} (km)');
    
    % Subplot 17: 5% Depth
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, Depth_data(i,2), color, num2str(i), sci_flag);
    quantity(i, t) = Depth_data(i,2);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('Depth_{5%} (km)');
    
    % Subplot 18: 95% Depth
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, Depth_data(i,3), color, num2str(i), sci_flag);
    quantity(i, t) = Depth_data(i,3);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('Depth_{95%} (km)');
    
    % Subplot 19: Fractal Dimension
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, df, color, num2str(i), sci_flag);
    quantity(i, t) = df;
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('{\itd_f}');
    
    % Subplot 20: Maximum Length
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, log10(Lmax), color, num2str(i), sci_flag);
    quantity(i, t) = log10(Lmax);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('log_{10}{\itL}_{max} (m)');
    
    % Subplot 21: Median Length
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, log10(Lmedian), color, num2str(i), sci_flag);
    quantity(i, t) = log10(Lmedian);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('log_{10}{\itL}_{median} (m)');
    
    % Subplot 22: Minimum Length
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, log10(Lmin), color, num2str(i), sci_flag);
    quantity(i, t) = log10(Lmin);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('log_{10}{\itL}_{min} (m)');
    
%     % Subplot 23: Volume
%     t = t + 1;
%     subplot(Numrow, Numcol, t); hold on;
%     plt_scatter(P_LR, log10(Vellipsoid), color, num2str(i), sci_flag);
%     quantity(i, t) = log10(Lmin);
%     set(gca, 'fontsize', 12);
%     xlabel('{\itp}_{LRT}'); ylabel('log_{10}{\itV}_{ell} (m^3)');
%     
    
    % Subplot 23: Volume
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, (Lmax)/ (Lmin), color, num2str(i), sci_flag);
    quantity(i, t) = (Lmax)/ (Lmin);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('{\itL}_{max}/{\itL}_{min}');
    
    % Subplot 23: Volume
    t = t + 1;
    subplot(Numrow, Numcol, t); hold on;
    plt_scatter(P_LR, log10(Vellipsoid), color, num2str(i), sci_flag);
    quantity(i, t) = log10(Vellipsoid);
    set(gca, 'fontsize', 12);
    xlabel('{\itp}_{LRT}'); ylabel('log_{10}{\itV}_{ell} (m^3)');  
    
    
end

% Store ylabels for later use
ylabelCell = cell(Numrow*Numcol, 1);
for i = 1:t
    subplot(Numrow, Numcol, i);
    ax = gca;
    ylabelCell{i} = ax.YLabel.String;
end

%% Create Figure 3 with nexttile
fig3 = figure('units', 'normalized', 'position', [0.05, 0.05, 0.4, 0.8]);
tiledlayout(Numrow, Numcol, 'TileSpacing', 'compact', 'Padding', 'compact');

P_LR = Case_Summary(:, 12);
TGR_5 = find(P_LR <= 0.05);
GR_5 = find(P_LR > 0.05);
TGR_10 = find(P_LR <= 0.1);
GR_10 = find(P_LR > 0.1);

len_TGR_5 = length(TGR_5);
len_TGR_10 = length(TGR_10);
num_rand = 1e3;

for j = 1:t
    nexttile;
    hold on;
    
    target_quantity = quantity(:, j);
    x_bins = unique(target_quantity);  
    
    % Random permutation tests
    mas_dis_5_rand = zeros(num_rand, 1);
    mas_dis_10_rand = zeros(num_rand, 1);
    cvm_dis_5_rand = zeros(num_rand, 1);
    cvm_dis_10_rand = zeros(num_rand, 1);

    for kk = 1:num_rand
        rand_indices = randperm(length(target_quantity));
        sel5 = rand_indices(1:len_TGR_5);
        rem5 = rand_indices(len_TGR_5+1:end);
        sel10 = rand_indices(1:len_TGR_10);
        rem10 = rand_indices(len_TGR_10+1:end);

        F1 = manual_ecdf(target_quantity(sel5), x_bins);
        F2 = manual_ecdf(target_quantity(rem5), x_bins);
        F3 = manual_ecdf(target_quantity(sel10), x_bins);
        F4 = manual_ecdf(target_quantity(rem10), x_bins);

        mas_dis_5_rand(kk) = max(abs(F1 - F2));
        mas_dis_10_rand(kk) = max(abs(F3 - F4));
        cvm_dis_5_rand(kk) = cramer_von_mises_stat(F1, F2);
        cvm_dis_10_rand(kk) = cramer_von_mises_stat(F3, F4);
    end

    % Real group comparison
    F_TGR5 = manual_ecdf(target_quantity(TGR_5), x_bins);
    F_GR5 = manual_ecdf(target_quantity(GR_5), x_bins);
    F_TGR10 = manual_ecdf(target_quantity(TGR_10), x_bins);
    F_GR10 = manual_ecdf(target_quantity(GR_10), x_bins);

    mas_dis_5 = max(abs(F_TGR5 - F_GR5));
    mas_dis_10 = max(abs(F_TGR10 - F_GR10));
    P_KS_5 = mean(mas_dis_5_rand >= mas_dis_5);
    P_KS_10 = mean(mas_dis_10_rand >= mas_dis_10);
    
    W2_5 = cramer_von_mises_stat(F_TGR5, F_GR5);
    W2_10 = cramer_von_mises_stat(F_TGR10, F_GR10);
    P_CvM_5 = mean(cvm_dis_5_rand >= W2_5);
    P_CvM_10 = mean(cvm_dis_10_rand >= W2_10);

    % Plot CDFs
    plot(F_TGR5, x_bins, 'r-', 'LineWidth', 1.5);
    plot(F_GR5, x_bins, 'b-', 'LineWidth', 1.5);
    plot(F_TGR10, x_bins, 'r--', 'LineWidth', 1.5);
    plot(F_GR10, x_bins, 'b--', 'LineWidth', 1.5);

    ylabel(ylabelCell{j});
    xlabel('CDF');
    set(gca, 'fontsize', 12);
    box off;
    
    % Improved title formatting
    title(sprintf('KS_{5%%}: %.2f  CvM_{5%%}: %.2f\nKS_{10%%}: %.2f  CvM_{10%%}: %.2f', ...
        P_KS_5, P_CvM_5, P_KS_10, P_CvM_10), 'FontSize', 10);
    
end

% After the for-loop ends
nexttile;
% Plot invisible lines for legend representation
h1 = plot(nan, nan, 'r-', 'LineWidth', 1.5); % TGR, p <= 0.05
hold on;
h2 = plot(nan, nan, 'b-', 'LineWidth', 1.5); % GR, p > 0.05
h3 = plot(nan, nan, 'r--', 'LineWidth', 1.5); % TGR, p <= 0.10
h4 = plot(nan, nan, 'b--', 'LineWidth', 1.5); % GR, p > 0.10

legend([h1, h2, h3, h4], ...
    {'TGR ({\itp}_{LRT} \leq 0.05)', 'GR ({\itp}_{LRT} > 0.05)', ...
     'TGR ({\itp}_{LRT} \leq 0.10)', 'GR ({\itp}_{LRT} > 0.10)'}, ...
     'Location', 'east', 'FontSize', 10);
axis off; % Hide the axis, only show legend
%}
%% Create Figure 4 (Excluding Experiments) with nexttile
fig4 = figure('units', 'normalized', 'position', [0.05, 0.05, 0.4, 0.8]);
tiledlayout(Numrow, Numcol, 'TileSpacing', 'compact', 'Padding', 'compact');

P_LR = Case_Summary(:, 12);
Sci_flags = Data_Flag(:,9);
jkf = find(Sci_flags == 0);
P_LR = P_LR(jkf);
quantity_filtered = quantity(jkf, :);

TGR_5 = find(P_LR <= 0.05);
GR_5 = find(P_LR > 0.05);
TGR_10 = find(P_LR <= 0.1);
GR_10 = find(P_LR > 0.1);

len_TGR_5 = length(TGR_5);
len_TGR_10 = length(TGR_10);

for j = 1:t
    nexttile;
    hold on;
    
    target_quantity = quantity_filtered(:, j);
    x_bins = unique(target_quantity);
    
    % Random permutation tests
    mas_dis_5_rand = zeros(num_rand, 1);
    mas_dis_10_rand = zeros(num_rand, 1);
    cvm_dis_5_rand = zeros(num_rand, 1);
    cvm_dis_10_rand = zeros(num_rand, 1);

    for kk = 1:num_rand
        rand_indices = randperm(length(target_quantity));
        sel5 = rand_indices(1:len_TGR_5);
        rem5 = rand_indices(len_TGR_5+1:end);
        sel10 = rand_indices(1:len_TGR_10);
        rem10 = rand_indices(len_TGR_10+1:end);

        F1 = manual_ecdf(target_quantity(sel5), x_bins);
        F2 = manual_ecdf(target_quantity(rem5), x_bins);
        F3 = manual_ecdf(target_quantity(sel10), x_bins);
        F4 = manual_ecdf(target_quantity(rem10), x_bins);

        mas_dis_5_rand(kk) = max(abs(F1 - F2));
        mas_dis_10_rand(kk) = max(abs(F3 - F4));
        cvm_dis_5_rand(kk) = cramer_von_mises_stat(F1, F2);
        cvm_dis_10_rand(kk) = cramer_von_mises_stat(F3, F4);
    end

    % Real group comparison
    F_TGR5 = manual_ecdf(target_quantity(TGR_5), x_bins);
    F_GR5 = manual_ecdf(target_quantity(GR_5), x_bins);
    F_TGR10 = manual_ecdf(target_quantity(TGR_10), x_bins);
    F_GR10 = manual_ecdf(target_quantity(GR_10), x_bins);

    mas_dis_5 = max(abs(F_TGR5 - F_GR5));
    mas_dis_10 = max(abs(F_TGR10 - F_GR10));
    P_KS_5 = mean(mas_dis_5_rand >= mas_dis_5);
    P_KS_10 = mean(mas_dis_10_rand >= mas_dis_10);
    
    W2_5 = cramer_von_mises_stat(F_TGR5, F_GR5);
    W2_10 = cramer_von_mises_stat(F_TGR10, F_GR10);
    P_CvM_5 = mean(cvm_dis_5_rand >= W2_5);
    P_CvM_10 = mean(cvm_dis_10_rand >= W2_10);

    % Plot CDFs
    plot(F_TGR5, x_bins, 'r-', 'LineWidth', 1.5);
    plot(F_GR5, x_bins, 'b-', 'LineWidth', 1.5);
    plot(F_TGR10, x_bins, 'r--', 'LineWidth', 1.5);
    plot(F_GR10, x_bins, 'b--', 'LineWidth', 1.5);

    ylabel(ylabelCell{j});
    xlabel('CDF');
    set(gca, 'fontsize', 12);
    box off;
    
    title(sprintf('KS_{5%%}: %.2f  CvM_{5%%}: %.2f\nKS_{10%%}: %.2f  CvM_{10%%}: %.2f', ...
        P_KS_5, P_CvM_5, P_KS_10, P_CvM_10), 'FontSize', 10);

end

% After the for-loop ends
nexttile;
% Plot invisible lines for legend representation
h1 = plot(nan, nan, 'r-', 'LineWidth', 1.5); % TGR, p <= 0.05
hold on;
h2 = plot(nan, nan, 'b-', 'LineWidth', 1.5); % GR, p > 0.05
h3 = plot(nan, nan, 'r--', 'LineWidth', 1.5); % TGR, p <= 0.10
h4 = plot(nan, nan, 'b--', 'LineWidth', 1.5); % GR, p > 0.10

legend([h1, h2, h3, h4], ...
    {'TGR ({\itp}_{LRT} \leq 0.05)', 'GR ({\itp}_{LRT} > 0.05)', ...
     'TGR ({\itp}_{LRT} \leq 0.10)', 'GR ({\itp}_{LRT} > 0.10)'}, ...
     'Location', 'east', 'FontSize', 10);
axis off; % Hide the axis, only show legend

% Helper functions
function [] = plt_scatter(X, Y, color, index, sci_flag)
    scatter(X, Y, 150, color, 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
    text(X, Y, index, 'Color', 'w', 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    if sci_flag == 1
        scatter(X, Y, 180, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
end

function custom_cmap = get_colormap()
    n1 = round(0.05 / 0.01);
    n2 = round(0.05 / 0.01);
    n3 = 1 + round(0.9 / 0.01);
    red = [0.85, 0.1, 0.1];
    orange = [1.0, 0.5, 0.0];
    yellow = [1.0, 0.9, 0.3];
    blue = [0.2, 0.4, 0.9];
    purple = [0.5, 0.2, 0.6];
    cmap1 = [linspace(red(1), orange(1), n1)', ...
             linspace(red(2), orange(2), n1)', ...
             linspace(red(3), orange(3), n1)'];
    cmap2 = [linspace(orange(1), yellow(1), n2)', ...
             linspace(orange(2), yellow(2), n2)', ...
             linspace(orange(3), yellow(3), n2)'];
    cmap3 = [linspace(blue(1), purple(1), n3)', ...
             linspace(blue(2), purple(2), n3)', ...
             linspace(blue(3), purple(3), n3)'];
    custom_cmap = [cmap1; cmap2; cmap3];
end

function F = manual_ecdf(data, x_bins)
    data = data(~isnan(data));
    n = length(data);
    if n == 0
        F = zeros(size(x_bins));
    else
        F = arrayfun(@(x) sum(data <= x) / n, x_bins);
    end
end

function W2 = cramer_von_mises_stat(F1, F2)
    W2 = sum((F1 - F2).^2);
end
