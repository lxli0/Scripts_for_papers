clc; clear; close all

%% ——— Figure (Nature-style footprint) ———
figure('Units','normalized','Position',[0.06 0.06 0.3 0.55]);

%% ——— Load file names ———
namelist  = dir('./Data/*.txt');
file_name = {namelist.name};
case_num  = numel(file_name);
sites     = cell(case_num,1);
for i = 1:case_num
    [~, sites{i}] = fileparts(file_name{i});
end

%% ——— Load summary and custom colormap ———
load('./Case_Summary.mat');
custom_cmap = get_colormap();
Data_Flag   = load('/Users/linxuanli/Desktop/IIS-GRvsTGR/Final/Observation/Data/_Data_Flag.md');
Depth_data  = load('/Users/linxuanli/Desktop/IIS-GRvsTGR/Final/Observation/Data/_Depth.md');

[~, sort_idx] = sort(Case_Summary(:,12),'descend');
sites         = sites(sort_idx);
Case_Summary  = Case_Summary(sort_idx,:);
Data_Flag     = Data_Flag(sort_idx,:);
Depth_data    = Depth_data(sort_idx,:);

%% ——— Prealloc ———
Obs_max_all   = nan(case_num,1);
logLmax_all   = nan(case_num,1);
logLmin_all   = nan(case_num,1);
logVellip_all = nan(case_num,1);
Depth_5_all   = nan(case_num,1);
Depth_median_all= nan(case_num,1);
Depth_95_all  = nan(case_num,1);
df_all        = nan(case_num,1);
P_LR          = nan(case_num,1);
r_ratio       = nan(case_num,1);
aspect_ratio  = nan(case_num,1);

%% ——— Loop: compute geometry metrics ———
for i = 1:case_num
    site = sites{i};
    D = load(['/Users/linxuanli/Desktop/IIS-GRvsTGR/Final/Observation/Data/' site '.txt']);
    D = sortrows(D,1);

    mc       = Case_Summary(i,1);
    Num      = Case_Summary(i,2);
    b_GR     = Case_Summary(i,3);
    b_TGR    = Case_Summary(i,6);
    m_corner = Case_Summary(i,9);
    P_LR(i)  = Case_Summary(i,12);
    Obs_max  = Case_Summary(i,15);
    a_GR     = log10(Num) + b_GR * mc; %#ok<NASGU>

    % threshold on Mc
    D  = D(D(:,5) >= mc - 1e-12, :);
    T0 = D(:,1);
    X0 = D(:,2) - mean(D(:,2));
    Y0 = D(:,3) - mean(D(:,3));
    Z0 = -abs(D(:,4) - mean(D(:,4)));
    M0 = D(:,5); %#ok<NASGU>

    flags = Data_Flag(i,:);  % [.. 4:deg 5:depth 6:time .. 9:scientific_flag]
    if flags(4) == 1, X0 = X0*110e3; Y0 = Y0*110e3; end   % deg→m
    if flags(5) == 1, Z0 = Z0*1e3; end                   % km→m
    if flags(6) == 2, T0 = (T0 - min(T0))*365.25; end    % yr→days (unused)

    if ~any(isnan([X0(:);Y0(:);Z0(:)]))
        [Lmax, Lmedian, Lmin, Vellip, ~] = get_spatial_ellipsoid(X0,Y0,Z0,2);
        df_all(i) = fractdime(X0,Y0,Z0,0);
    else
        Lmax = nan; Lmin = nan; Vellip = nan; df_all(i) = nan;
    end

    % store
    Obs_max_all(i)   = Obs_max;
    logLmax_all(i)   = log10(Lmax);
    logLmin_all(i)   = log10(Lmin);
    logVellip_all(i) = log10(Vellip);
    Depth_5_all(i)   = Depth_data(i,2);
    Depth_median_all(i)   = Depth_data(i,1);
    Depth_95_all(i)  = Depth_data(i,3);

    % log2 aspect ratio and log10 size ratio
    aspect_ratio(i) = log2(Lmax/Lmin);
    r_1MPa    = (7/16/1e6 * 10^(1.5*(Obs_max+6.07))).^(1/3);
    r_ratio(i)= log10(r_1MPa / Lmax);
end

%% ——— Grouping (TGR vs GR) ———
Sci_flags = Data_Flag(:,9);

% For most panels keep the "good-quality" subset
good_idx  = (Sci_flags == 0);

P_LR_f        = P_LR(good_idx);
Depth5_f      = Depth_5_all(good_idx);
Depth95_f     = Depth_95_all(good_idx);
Depth_median_f = Depth_median_all(good_idx);
Lmax_f        = logLmax_all(good_idx);
Vellip_f      = logVellip_all(good_idx);
rratio_f      = r_ratio(good_idx);
aspect_ratio_f= aspect_ratio(good_idx);

isTGR_f   = P_LR_f <= 0.05;
G_f       = [repmat({'TGR'},sum(isTGR_f),1); repmat({'GR'},sum(~isTGR_f),1)];
COL       = [0.85, 0.1, 0.1;0.2, 0.4, 0.9];

% For the df panel: include ALL cases
valid_df  = ~isnan(df_all) & ~isnan(P_LR);
df_allv   = df_all(valid_df);
P_LR_allv = P_LR(valid_df);
isTGR_df  = P_LR_allv <= 0.05;
G_df      = [repmat({'TGR'},sum(isTGR_df),1); repmat({'GR'},sum(~isTGR_df),1)];

%% ——— (1) Seismicity depths (start & end, overlaid) ———
ax = subplot(2,3,1);
y  = [Depth_median_f(isTGR_f); Depth_median_f(~isTGR_f)];
plot_median_iqr(y, G_f, COL);

ylabel('Seismicity depth (km)');
ylim([2 5]);
set(ax, 'YDir', 'reverse');    % Reverse y-axis direction
set(ax, 'XAxisLocation', 'top'); % Move x-axis to top
hold off
nature_axes(ax);

%% ——— (2) Seismic cloud length (log10 scale values, pseudo 10^ labels) ———
ax = subplot(2,3,2);
y  = [Lmax_f(isTGR_f); Lmax_f(~isTGR_f)]-3;  % log10(km)
plot_median_iqr(y, G_f, COL);
ylabel('Seismic cloud length (km)');
ylim([-0.5 1.5]);
% integer pseudo 10^ labels
yticks(0:1:1); 
yticklabels(arrayfun(@(t) sprintf('10^{%d}', round(t)), yticks, 'UniformOutput', false));
nature_axes(ax);

%% ——— (3) Seismic cloud volume (converted to km^3, pseudo 10^ labels) ———
ax = subplot(2,3,3);
y  = [Vellip_f(isTGR_f); Vellip_f(~isTGR_f)] - 9;  % convert m³→km³
plot_median_iqr(y, G_f, COL);
ylabel('Seismic cloud volume (km^3)');
ylim([-0.5 3.5]);
yticks(-0:1:3);
yticklabels(arrayfun(@(t) sprintf('10^{%d}', round(t)), yticks, 'UniformOutput', false));
nature_axes(ax);


%% ——— (4) Aspect ratio (log2 values, pseudo 2^ labels) ———
ax = subplot(2,3,4);
y  = [aspect_ratio_f(isTGR_f); aspect_ratio_f(~isTGR_f)];  % values are log2(Lmax/Lmin)
plot_median_iqr(y, G_f, COL);
ylabel({'Aspect ratio: cloud length','divided by width'});
% choose neat integer powers for 2^· labels
ylim([0.8 4.2]);
ticks = 1:1:4;
if numel(ticks) < 3
    ticks = round(yl(1)):round(yl(2)); % fallback
end
yticks(ticks);
yticklabels(arrayfun(@(t) sprintf('2^{%d}', round(t)), yticks, 'UniformOutput', false));
nature_axes(ax);


%% ——— (5) Relative size: cloud length / largest-event radius (log10 values) ———
ax = subplot(2,3,5);
plot([0 3],[0 0],'--','color',[0.5 0.5 0.5]);
hold on; 
y  = [rratio_f(isTGR_f); rratio_f(~isTGR_f)]; % log10(Lmax / r_1MPa)
plot_median_iqr(y, G_f, COL);
ylabel({'Largest-event radius divided','by seismic cloud length'});
% integer pseudo 10^ labels
ylim([-2 0.2]);
yticks([-2 -1 0]);
yticklabels(arrayfun(@(t) sprintf('10^{%d}', round(t)), yticks, 'UniformOutput', false));
nature_axes(ax);



%% ——— (6) Fractal dimension ———
ax = subplot(2,3,6);
y  = [df_allv(isTGR_df); df_allv(~isTGR_df)];
plot_median_iqr(y, G_df, COL);
ylabel({'Fractal dimension', 'of pairwise distance'});
ylim([1 3]);
nature_axes(ax);

%% ——— Helpers ———
function nature_axes(ax)
% Apply a “Nature-style” axis treatment.
    set(ax, 'TickDir','out', ...
            'TickLength',[0.015 0.015], ...
            'LineWidth',1, ...
            'Box','off', ...
            'FontName','Helvetica', ...
            'FontSize',13, ...
            'XMinorTick','off', ...
            'YMinorTick','on');
end

function custom_cmap = get_colormap()
    n1 = round(0.05/0.01);
    n2 = round(0.05/0.01);
    n3 = 1 + round(0.9/0.01);
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

function plot_median_iqr(y, G, COL, cap_length)
% Plot median (dot) and IQR (vertical line with small caps) for each group in G.
% Optional: cap_length controls the horizontal cap length (0.1–0.3 recommended).
    if nargin < 4, cap_length = 0.24; end
    hold on
    % normalize labels
    if iscell(G) || ischar(G)
        gvals = unique(G,'stable');        sel = @(k) strcmp(G,gvals{k}); xtlabs = string(gvals);
    elseif isstring(G)
        gvals = unique(G,'stable');        sel = @(k) (G==gvals(k));      xtlabs = gvals;
    elseif iscategorical(G)
        gvals = categories(G);             sel = @(k) (G==categorical(gvals{k})); xtlabs = string(gvals);
    else
        gvals = unique(G,'stable');        sel = @(k) (G==gvals(k));      xtlabs = string(gvals);
    end
    ng = numel(gvals);
    for k = 1:ng
        ii=ng+1-k;
        idx = sel(k);
        d   = y(idx); d = d(~isnan(d));
        if isempty(d), continue; end
        p25 = prctile(d,25); med = median(d); p75 = prctile(d,75);
        c = COL(1+mod(k-1,size(COL,1)),:);
        % IQR + caps
        line([ii ii],[p25 p75],'Color',c,'LineWidth',2.5);
        line([ii-cap_length/2 ii+cap_length/2],[p25 p25],'Color',c,'LineWidth',2.5);
        line([ii-cap_length/2 ii+cap_length/2],[p75 p75],'Color',c,'LineWidth',2.5);
        % median
        plot(ii, med, 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'MarkerSize', 12);
    end
    xlim([0.7 ng+0.3]); xticks(1:ng); xticklabels(flip(xtlabs));
    box off; hold off
end

