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

% === Load case summary and data flags ===
load('Case_Summary.mat');
Data_Flag = load('./Data/_Data_Flag.md');   % assumes numeric table in md
[~, sort_idx] = sort(Case_Summary(:,12), 'descend');
sites = sites(sort_idx);
Case_Summary = Case_Summary(sort_idx, :);
Data_Flag = Data_Flag(sort_idx,:);

% === Backup for subset ===
Case_Summary_all = Case_Summary;
Data_Flag_all    = Data_Flag;

sci_flag      = Data_Flag_all(:,9);
filtered_idx  = find(sci_flag ~= 1);
Case_Summary  = Case_Summary_all(filtered_idx, :);

% === Create figure with 4x3 subplots (large) ===
% Histogram bins
bin_edges   = -1:1:6;
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% ---------- Data for filtered cases ----------
P_LR       = Case_Summary(:,12);
m_corner   = Case_Summary(:,9);
mmax_obs   = Case_Summary(:,15);
mode_max_GR= Case_Summary(:,16); 
mc         = Case_Summary(:,1);  
Num        = Case_Summary(:,2); 

TGR_5  = find(P_LR <= 0.05); GR_5  = find(P_LR > 0.05);
TGR_10 = find(P_LR <= 0.10); GR_10 = find(P_LR > 0.10);

fig1 = figure('units', 'normalized', 'position', [0.05, 0.05, 0.85, 0.9]);

% --- Subplot 9: m_corner (filtered) ---
hmcor = subplot(3, 3, 1);
[counts1_TGR_5,  ~] = histcounts(m_corner(TGR_5),  bin_edges);
[counts1_GR_5,   ~] = histcounts(m_corner(GR_5),   bin_edges);
[counts1_TGR_10, ~] = histcounts(m_corner(TGR_10), bin_edges);
[counts1_GR_10,  ~] = histcounts(m_corner(GR_10),  bin_edges);

plot(bin_centers, counts1_TGR_5,  '-o', 'Color', [0.85, 0.10, 0.10], 'LineWidth', 2); hold on;
plot(bin_centers, counts1_GR_5,   '-o', 'Color', [0.20, 0.40, 0.90], 'LineWidth', 2);
plot(bin_centers, counts1_TGR_10, '--o', 'Color', [0.85, 0.10, 0.10], 'LineWidth', 2);
plot(bin_centers, counts1_GR_10,  '--o', 'Color', [0.20, 0.40, 0.90], 'LineWidth', 2);
ylim([0 8]);
xlim([-1 6]);
xlabel('Estimated {\itm}_{corner}','Interpreter','tex'); ylabel('Count');
set(gca, 'FontSize', 14, 'Box','off');
pos = get(hmcor, 'Position');
pos(4) = pos(4) - 0.01; % widen
set(hmcor, 'Position', pos);

% --- Subplot 12: m_max (filtered) ---
hmax = subplot(3, 3, 2);
[counts2_TGR_5,  ~] = histcounts(mmax_obs(TGR_5),  bin_edges);
[counts2_GR_5,   ~] = histcounts(mmax_obs(GR_5),   bin_edges);
[counts2_TGR_10, ~] = histcounts(mmax_obs(TGR_10), bin_edges);
[counts2_GR_10,  ~] = histcounts(mmax_obs(GR_10),  bin_edges);

plot(bin_centers, counts2_TGR_5,  '-o', 'Color', [0.85, 0.10, 0.10], 'LineWidth', 2); hold on;
plot(bin_centers, counts2_GR_5,   '-o', 'Color', [0.20, 0.40, 0.90], 'LineWidth', 2);
plot(bin_centers, counts2_TGR_10, '--o', 'Color', [0.85, 0.10, 0.10], 'LineWidth', 2);
plot(bin_centers, counts2_GR_10,  '--o', 'Color', [0.20, 0.40, 0.90], 'LineWidth', 2);
ylim([0 8]);
xlim([-1 6]);
xlabel('Observed {\itm}_{max}','Interpreter','tex'); ylabel('Count');
set(gca, 'FontSize', 14, 'Box','off');
pos = get(hmax, 'Position');
pos(4) = pos(4) - 0.01; % widen
set(hmax, 'Position', pos);

% =================== POLISHING & SHARED LEGEND ===================

% Nature-ish axis styling (keeps your large layout)
set([hmcor hmax], 'FontName','Arial', 'FontSize', 14, 'LineWidth',0.75, ...
    'TickDir','out', 'Box','off');

% Build ONE shared legend using the four lines from hmcor
axChildren = findobj(hmcor,'Type','line');   % returns in reverse order
axChildren = flipud(axChildren);             % restore plotting order
% Safety: ensure we have at least 4 plotted lines
if numel(axChildren) >= 4
    h1 = axChildren(1); h2 = axChildren(2); h3 = axChildren(3); h4 = axChildren(4);
else
    error('Expected four line objects in hmcor for legend creation.');
end

hL = legend(hmcor, [h1 h2 h3 h4], ...
   'TGR ({\itp}_{LRT} \leq 0.05)', ...
   'GR ({\itp}_{LRT} > 0.05)', ...
   'TGR ({\itp}_{LRT} \leq 0.10)', ...
   'GR ({\itp}_{LRT} > 0.10)', ...
   'Interpreter','tex', 'NumColumns',4, 'Box','off');
set(hL, 'Units','normalized');

% Place legend just outside the two panels (right side)
p1 = get(hmcor,'Position');   % [x y w h] normalized
p2 = get(hmax, 'Position');

bbox_left   = min(p1(1), p2(1));
bbox_right  = max(p1(1)+p1(3), p2(1)+p2(3));
bbox_bottom = min(p1(2), p2(2));
bbox_top    = max(p1(2)+p1(4), p2(2)+p2(4));
bbox_height = bbox_top - bbox_bottom;

legW = 0.18;                        % legend width (tweak if needed)
legH = 0.15;  % height cap for neatness
gap  = 0.02;                        % gap between right axis and legend

legLeft = bbox_left+0.15;
legBot  = bbox_top-0.05;%bbox_bottom + (bbox_height - legH)/2; % vertical center
set(hL, 'Position', [legLeft, legBot, legW, legH]);

%
% ===== Nature-style single panel (histogram curves as probabilities) =====
% Fallbacks if edges/centers aren't defined in your workspace:
% if ~exist('bin_edges','var'),   bin_edges = -1:1:6; end
% if ~exist('bin_centers','var'), bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2; end

% --- Figure + defaults ---
fig = figure('Units','normalized','Position',[0.05, 0.05, 0.25, 0.2]);

set(fig, 'DefaultAxesFontName','Arial', ...
         'DefaultAxesFontSize',14, ...
         'DefaultAxesLineWidth',0.75, ...
         'DefaultLineLineWidth',1.25, ...
         'DefaultLineMarkerSize',4);

ax = axes(fig); hold(ax,'on');
set(ax,'Box','off','TickDir','out');

% --- Compute histograms ---
[counts1_aoverbminus1,  ~] = histcounts(mode_max_GR - 1.0, bin_edges);
[counts1_aoverbminus05, ~] = histcounts(mode_max_GR - 0.5, bin_edges);

% Safe normalization to probabilities
safe_norm = @(v) (sum(v)>0) .* (v / max(sum(v), eps));

p_aoverbminus05 = safe_norm(counts1_aoverbminus05);
p_aoverbminus1  = safe_norm(counts1_aoverbminus1);
p_TGR5          = safe_norm(counts1_TGR_5);

% --- Color/line scheme (colorblind/print-friendly) ---
% --- Color/line scheme (colorblind- and print-friendly warm palette) ---
col_orange = [0.90 0.50 0.00];   % bright orange
col_purple = [0.55 0.30 0.65];   % muted purple (contrasts both warm tones)
col_red    = [0.85, 0.10, 0.10];  % strong red

hold on;
% --- Plot (Nature-style) ---
h1 = plot(ax, bin_centers, p_aoverbminus05, '-o', ...
          'Color', col_orange, 'MarkerFaceColor','none', ...
          'MarkerEdgeColor', col_orange, 'LineWidth',2.0);

h2 = plot(ax, bin_centers, p_aoverbminus1,  '-o', ...
          'Color', col_purple, 'MarkerFaceColor','none', ...
          'MarkerEdgeColor', col_purple, 'LineWidth',2.0);

h3 = plot(ax, bin_centers, p_TGR5, '-o', ...
          'Color', col_red, 'MarkerFaceColor','none', ...
          'MarkerEdgeColor', col_red, 'LineWidth',2.0);

      % --- Labels/limits ---
xlabel(ax, 'Value','Interpreter','tex');
ylabel(ax, 'Frequency','Interpreter','tex');
xlim(ax, [bin_edges(1) bin_edges(end)]);

% Optional minor grid for readability (very light)
% set(ax,'YMinorTick','on');

% --- Legend outside to the right ---
lgd = legend(ax, [h1 h2 h3], ...
    '{\ita}_{GR}/{\itb}_{GR} - 0.5', ...
    '{\ita}_{GR}/{\itb}_{GR} - 1', ...
    'Estimated {\itm}_{corner}^{TGR}', ...
    'Interpreter','tex','NumColumns',1,'Box','off','Location','northeast');
ylim([0 0.6]);
% Optional panel letter
% text(ax, 0.0, 1.05, 'c', 'Units','normalized','FontWeight','bold');

