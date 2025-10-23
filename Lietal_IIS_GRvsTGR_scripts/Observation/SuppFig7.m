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

[~, sort_idx] = sort(Case_Summary(:, 12), 'descend');
sites = sites(sort_idx);
Case_Summary = Case_Summary(sort_idx, :);
Data_Flag = Data_Flag(sort_idx, :);
Depth_data = Depth_data(sort_idx, :);

Sci_flags = Data_Flag(:, 9);
jkf = find(Sci_flags == 0);
sites = sites(jkf);
Case_Summary = Case_Summary(jkf, :);
Data_Flag = Data_Flag(jkf, :);
Depth_data = Depth_data(jkf, :);
case_num = length(sites);

% === Create figures ===
fig1 = figure('units', 'normalized', 'position', [0.05, 0.05, 0.4, 0.4]);
fig2 = figure('units', 'normalized', 'position', [0.05, 0.05, 0.5, 0.4]);

Numrow1 = 2; Numcol1 = 3;
Numrow2 = 2; Numcol2 = 4;

a_GRs = zeros(case_num, 1);
Obs_max_all = zeros(case_num, 1);
m_corner_all = zeros(case_num, 1);
logLmax_all = zeros(case_num, 1);
logLmedian_all = zeros(case_num, 1);
logLmin_all = zeros(case_num, 1);
logVellip_all = zeros(case_num, 1);

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
    a_GR = log10(Num) + b_GR * mc;

    D = D(D(:,5) >= mc - 1e-12, :);
    T0 = D(:,1);
    X0 = D(:,2) - mean(D(:,2));
    Y0 = D(:,3) - mean(D(:,3));
    Z0 = -abs(D(:,4) - mean(D(:,4)));
    M0 = D(:,5);

    flags = Data_Flag(i,:);
    if flags(4) == 1, X0 = X0 * 110e3; Y0 = Y0 * 110e3; end
    if flags(5) == 1, Z0 = Z0 * 1e3; end
    if flags(6) == 2, T0 = (T0 - min(T0)) * 365.25; end

    if ~isnan(X0)
        [Lmax, Lmedian, Lmin, Vellipsoid, ~] = get_spatial_ellipsoid(X0, Y0, Z0, 2);
    else
        Lmax = nan; Lmedian = nan; Lmin = nan; Vellipsoid = nan;
    end

    a_GRs(i) = a_GR;
    Obs_max_all(i) = Obs_max;
    m_corner_all(i) = m_corner;
    logLmax_all(i) = log10(Lmax);
    logLmedian_all(i) = log10(Lmedian);
    logLmin_all(i) = log10(Lmin);
    logVellip_all(i) = log10(Vellipsoid);

    color = custom_cmap(1 + floor(round(P_LR, 2) * 100), :);
    sci_flag = flags(9);

    figure(fig1);
    t = 1; subplot(Numrow1, Numcol1, t); hold on;
    plt_scatter(a_GR, Obs_max, color, num2str(i), sci_flag);

    t = 2; subplot(Numrow1, Numcol1, t); hold on;
    plt_scatter(a_GR, m_corner, color, num2str(i), sci_flag);

    t = 3; subplot(Numrow1, Numcol1, t); hold on;
    plt_scatter(a_GR, log10(Lmax), color, num2str(i), sci_flag);

    t = 4; subplot(Numrow1, Numcol1, t); hold on;
    plt_scatter(a_GR, log10(Lmedian), color, num2str(i), sci_flag);

    t = 5; subplot(Numrow1, Numcol1, t); hold on;
    plt_scatter(a_GR, log10(Lmin), color, num2str(i), sci_flag);

    t = 6; subplot(Numrow1, Numcol1, t); hold on;
    plt_scatter(a_GR, log10(Vellipsoid), color, num2str(i), sci_flag);

    figure(fig2);
    t = 1; subplot(Numrow2, Numcol2, t); hold on;
    plt_scatter(log10(Lmax), Obs_max, color, num2str(i), sci_flag);
    xlabel('log_{10} {\itL}_{max} (m)'); ylabel('Observed {\itm}_{max}'); ylim([-1 6]);

    t = 2; subplot(Numrow2, Numcol2, t); hold on;
    plt_scatter(log10(Lmedian), Obs_max, color, num2str(i), sci_flag);
    xlabel('log_{10} {\itL}_{median} (m)'); ylabel('Observed {\itm}_{max}'); ylim([-1 6]);

    t = 3; subplot(Numrow2, Numcol2, t); hold on;
    plt_scatter(log10(Lmin), Obs_max, color, num2str(i), sci_flag);
    xlabel('log_{10} {\itL}_{min} (m)'); ylabel('Observed {\itm}_{max}'); ylim([-1 6]);

    t = 4; subplot(Numrow2, Numcol2, t); hold on;
    plt_scatter(log10(Vellipsoid), Obs_max, color, num2str(i), sci_flag);
    xlabel('log_{10} {\itV}_{ellipsoid} (m^3)'); ylabel('Observed {\itm}_{max}'); ylim([-1 6]);

    t = 5; subplot(Numrow2, Numcol2, t); hold on;
    plt_scatter(log10(Lmax), m_corner, color, num2str(i), sci_flag);
    xlabel('log_{10} {\itL}_{max} (m)'); ylabel('{\itm}_{corner}'); ylim([-1 6]);

    t = 6; subplot(Numrow2, Numcol2, t); hold on;
    plt_scatter(log10(Lmedian), m_corner, color, num2str(i), sci_flag);
    xlabel('log_{10} {\itL}_{median} (m)'); ylabel('{\itm}_{corner}'); ylim([-1 6]);

    t = 7; subplot(Numrow2, Numcol2, t); hold on;
    plt_scatter(log10(Lmin), m_corner, color, num2str(i), sci_flag);
    xlabel('log_{10} {\itL}_{min} (m)'); ylabel('{\itm}_{corner}'); ylim([-1 6]);

    t = 8; subplot(Numrow2, Numcol2, t); hold on;
    plt_scatter(log10(Vellipsoid), m_corner, color, num2str(i), sci_flag);
    xlabel('log_{10} {\itV}_{ellipsoid} (m^3)'); ylabel('{\itm}_{corner}'); ylim([-1 6]);
end

% === Optional: Perform linear fitting after loop ===
figure(fig1);
Y_all = {Obs_max_all, m_corner_all, logLmax_all, logLmedian_all, logLmin_all, logVellip_all};
for t = 1:6
    subplot(Numrow1, Numcol1, t);
    xlim([0 6]);
    [~, ~, R2] = linear_fit_plot(a_GRs, Y_all{t});
    title(sprintf('R^2 = %.2f', R2));
    switch t
        case 1
            xlabel('{\ita}_{GR}'); ylabel('Observed {\itm}_{max}');
        case 2
            xlabel('{\ita}_{GR}'); ylabel('{\itm}_{corner}');
        case 3
            xlabel('{\ita}_{GR}'); ylabel('log_{10} {\itL}_{max}');
        case 4
            xlabel('{\ita}_{GR}'); ylabel('log_{10} {\itL}_{median}');
        case 5
            xlabel('{\ita}_{GR}'); ylabel('log_{10} {\itL}_{min}');
        case 6
            xlabel('{\ita}_{GR}'); ylabel('log_{10} {\itV}_{ellipsoid}');
    end
end

figure(fig2);
for i = 1:8
    subplot(Numrow2, Numcol2, i);
    xlim([8 16])
    ylim([-1 6])
    if i ~= 4 && i ~= 8
        subplot(Numrow2, Numcol2, i);
        lgpltL = 0:1:8
        pltL = 10.^lgpltL;
        for factor = -1:2:7
            m_pred = 2/3 * log10(pltL.^3 * 10.^factor * 16/7) - 6.07;
            plot(lgpltL, m_pred, 'k', 'linewidth', 1);
        end
        xlim([2 6]);
    end
end

function [] = plt_scatter(X, Y, color, index, sci_flag)
    scatter(X, Y, 150, color, 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
    set(gca,'fontsize',14)
end

function custom_cmap = get_colormap()
    n1 = round(0.05 / 0.01);
    n2 = round(0.05 / 0.01);
    n3 = 1 + round(0.9 / 0.01);
    red     = [0.85, 0.1, 0.1];
    orange  = [1.0, 0.5, 0.0];
    yellow  = [1.0, 0.9, 0.3];
    blue    = [0.2, 0.4, 0.9];
    purple  = [0.5, 0.2, 0.6];
    cmap1 = [linspace(red(1), orange(1), n1)', linspace(red(2), orange(2), n1)', linspace(red(3), orange(3), n1)'];
    cmap2 = [linspace(orange(1), yellow(1), n2)', linspace(orange(2), yellow(2), n2)', linspace(orange(3), yellow(3), n2)'];
    cmap3 = [linspace(blue(1), purple(1), n3)', linspace(blue(2), purple(2), n3)', linspace(blue(3), purple(3), n3)'];
    custom_cmap = [cmap1; cmap2; cmap3];
end

function [xfit, yfit, R2] = linear_fit_plot(X, Y)
    valid = ~(isnan(X) | isnan(Y));
    X = X(valid); Y = Y(valid);
    if length(X) <= 1 || length(unique(X)) == 1
        xfit = []; yfit = []; R2 = NaN;
        return;
    end
    p = polyfit(X, Y, 1);
    xfit = linspace(min(X), max(X), 100);
    yfit = polyval(p, xfit);
    hold on; plot(xfit, yfit, 'k--', 'LineWidth', 2.5);
    Yfit = polyval(p, X);
    R2 = 1 - sum((Y - Yfit).^2) / sum((Y - mean(Y)).^2);
end


