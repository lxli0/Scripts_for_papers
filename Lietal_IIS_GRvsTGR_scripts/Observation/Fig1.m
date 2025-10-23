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


[~, sort_idx] = sort(Case_Summary(:,12), 'descend');
sites = sites(sort_idx);
Case_Summary = Case_Summary(sort_idx,:);

% === Create figure ===
fig1=figure('units', 'normalized', 'position', [0.05, 0.05, 0.75, 0.9]);
fig2=figure('units', 'normalized', 'position', [0.05, 0.05, 0.75, 0.9]);
% === GR mmax comparison ===
figure(fig1);

figure(fig2);
subplot(3, 3, 1); 
maglim_min = -4; maglim_max = 6;
plot([maglim_min, maglim_max], [0 0], '--', ...
    'linewidth', 2, 'color', [0.5, 0.5, 0.5]); hold on;
grid off; box on;
xlim([maglim_min, maglim_max]);
ylim([-2.5, 2.5]);
xticks(maglim_min:1:maglim_max);
yticks(-3:1:3);
ylabel('{\itm}_{\rm max}^{\rm GR} - {\itm}_{\rm max}^{Observed}');
xlabel('{\itm}_{\rm max}^{Observed}');
set(gca, 'FontSize', 16, 'FontName', 'Arial', 'Box', 'off', 'TickDir', 'out');

% === TGR mmax comparison ===
subplot(3, 3, 2); 
plot([maglim_min, maglim_max], [0 0], '--', ...
    'linewidth', 2, 'color', [0.5, 0.5, 0.5]); hold on;
grid off; box on;
xlim([maglim_min, maglim_max]);
ylim([-2.5, 2.5]);
xticks(maglim_min:1:maglim_max);
yticks(-3:1:3);
ylabel('{\itm}_{\rm max}^{\rm TGR} - {\itm}_{\rm max}^{Observed}');
xlabel('{\itm}_{\rm max}^{Observed}');
set(gca, 'FontSize', 16, 'FontName', 'Arial', 'Box', 'off', 'TickDir', 'out');

alpha_GR_all = [];
alpha_TGR_all = [];

% === Loop over cases ===
for i = 1:case_num
    site = sites{i};
    file = ['./Data/', site, '.txt'];
    D = load(file); 
    D = sortrows(D, 1);
    D = D(D(:,1) > -1e4, :);
    T0 = D(:,1);
    M0 = D(:,5);

    mc = Case_Summary(i,1);
    P_LR = Case_Summary(i,12);
    alpha_GR = Case_Summary(i,13);
    alpha_TGR = Case_Summary(i,14);
    Obs_max = Case_Summary(i,15);
    mode_max_GR = Case_Summary(i,16);
    max_GR_5 = Case_Summary(i,17);
    max_GR_95 = Case_Summary(i,18);
    mode_max_TGR = Case_Summary(i,19);
    max_TGR_5 = Case_Summary(i,20);
    max_TGR_95 = Case_Summary(i,21);
    color = custom_cmap(1 + floor(round(P_LR, 2) * 100), :);
    color_save(i,:)=color;

    figure(fig1);
    subplot(2, 3, [1,2,3]); 
    idx = M0 - mc >= -1e-11;
    Dm = M0(idx) - mc;
    diffm = diff(M0);
    nonzero_elements = diffm(diffm ~= 0);
    delta = min(abs(nonzero_elements)) / 2;
    pl(i)=plt_MFD(Dm, color, num2str(i));
    leg{i}=[num2str(i),'-',site];
    
    
    figure(fig2);
    subplot(3, 3, 1); hold on;
    plot([Obs_max, Obs_max], [max_GR_5-Obs_max, max_GR_95-Obs_max] , 'color', color, 'linewidth', 2);
    scatter(Obs_max, mode_max_GR-Obs_max,150, color, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    text(Obs_max,mode_max_GR-Obs_max, num2str(i), 'Color', 'w', 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

    % === TGR mmax ===
    subplot(3,3,2); hold on;
    plot([Obs_max, Obs_max], [max_TGR_5-Obs_max, max_TGR_95-Obs_max], 'color', color, 'linewidth', 2);
    scatter(Obs_max, mode_max_TGR-Obs_max,150, color, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    text(Obs_max, mode_max_TGR-Obs_max, num2str(i), 'Color', 'w', 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% === Final formatting for MFD ===
figure(fig1);
subplot(2, 3, [1,2,3]); 
grid off; box on;
set(gca, 'FontSize', 16, 'FontName', 'Arial', 'Box', 'off', 'TickDir', 'out');
xlabel('Magnitude relative to completeness, {\itm} - {\itm}_C');
ylabel('{\itP}(≥{\it m})');
xlim([0 4.5]);
ylim([6e-5 1]);
pos = get(gca, 'Position'); 
pos(2)=pos(2)-0.05;
pos(4)=pos(4)+0.05;
set(gca, 'Position', pos);


% === Colorbar ===
colormap(custom_cmap);
cb = colorbar('north');
cb.Label.String = '{\itp}_{\rm LRT}';
cb.Label.FontSize = 16;
cb.Ticks = [0,0.1,0.2,0.6,1.0];
cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), cb.Ticks, 'UniformOutput', false);
cb.FontSize = 16;
% Make the colorbar shorter
cb.Position(1) = cb.Position(1) + 0.15;  % Shift to the right
cb.Position(2) = cb.Position(2) + 0.02;  % Shift to the right
cb.Position(3) = cb.Position(3) - 0.4;  % Reduce width
% Move the label to the top
cb.Label.VerticalAlignment = 'top';  % So it appears above the bar
cb.Label.HorizontalAlignment = 'center';
cb.Label.Position(2) = cb.Label.Position(2) + 0.6;  % Adjust upward offset as needed
% === inset ===
% Get the position of the current subplot
mainPos = get(gca, 'Position');  % [left bottom width height]

% Define inset position relative to subplot (e.g., upper right corner)
insetRel = [0.78, 0.72, 0.23, 0.34];  % [x0 y0 width height] relative to subplot

% Convert to figure units
insetAbs = [ ...
    mainPos(1) + insetRel(1) * mainPos(3), ...
    mainPos(2) + insetRel(2) * mainPos(4), ...
    insetRel(3) * mainPos(3), ...
    insetRel(4) * mainPos(4)];

% Create inset axes
insetAx1 = axes('Position', insetAbs);
CI;

    
subplot(2,2,3)
inset_map(sites,color_save);
pos = get(gca, 'Position'); 
pos(1)=pos(1)-0.04;
pos(2)=pos(2)-0.05;
pos(3)=pos(3)+0.15;
set(gca, 'Position', pos);




% === Legend at subplot(2,3,3) position ===
subplot(2,3,6); 
 pos = get(gca, 'Position'); 
pos(1) = pos(1) -0.05; % shift right by 0.02 units
pos(2) = pos(2)-0.03; % shift right by 0.02 units
set(gca, 'Position', pos);

delete(gca);
lgd = legend(pl, leg, 'NumColumns', 2);
set(lgd, 'Box', 'off', 'Position', pos, ...
    'FontSize', 14, 'FontName', 'Helvetica', 'Interpreter', 'none');

figure(fig2);
%%
hGR = subplot(3, 3, 1);
pos = get(hGR, 'Position');
%pos(1) = pos(1) -0.01; % shift right by 0.02 units
set(hGR, 'Position', pos);

hold on;


text(-3.8, -1.5, {'Small-scale', 'experiments'},'FontSize', 16);

% legend
h_mode = plot(nan, nan, 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 0.5);
h_ci = plot(nan, nan, 'k', 'LineWidth', 2);

lgd=legend([h_mode, h_ci], {'mode', '90% CI'}, ...
    'Location', 'northwest', 'FontSize', 16, 'Box', 'off');
pos = get(lgd, 'Position');
pos(2) = pos(2) +0.05; % shift right by 0.02 units
set(lgd, 'Position', pos);


hTGR=subplot(3, 3, 2);
pos = get(hTGR, 'Position');
%pos(1) = pos(1) + 0.05; % shift right by 0.02 units
set(hTGR, 'Position', pos);

text(-3.8, -1.5, {'Small-scale', 'experiments'},'FontSize', 16);

subplot(3,3,3)
% pos = get(gca, 'Position'); 
% pos(4)=pos(4)+0.10;
% set(gca, 'Position', pos);
plt_mcorner_his;

%run('/Users/linxuanli/Desktop/IIS-GRvsTGR/Final/SeismicCloud_vs_mmax/plt_radius.m');

%%
% Create figure
fig3 = figure('Units', 'normalized', 'Position', [0.05, 0.05, 0.75, 0.9]);

% Create subplot (2,2,3)
subplot(2,2,3)
pos = get(gca, 'Position'); 
pos(1)=pos(1)-0.04;
pos(2)=pos(2)-0.05;
pos(3)=pos(3)+0.15;
set(gca, 'Position', pos);

% --- Hide everything about the axes ---
axis off;             % hides box, ticks, labels, and grid
box off;              % ensures frame is gone too

% If you still need to draw something in it, keep its handle:
mainAx = gca;
set(mainAx, 'Color', 'w');   % optional: white background

mainPos = get(mainAx, 'Position');  % for inset placement

% === First inset ===
insetRel = [0.15, 0.15, 0.25, 0.3];
insetAbs = [ ...
    mainPos(1) + insetRel(1) * mainPos(3), ...
    mainPos(2) + insetRel(2) * mainPos(4), ...
    insetRel(3) * mainPos(3), ...
    insetRel(4) * mainPos(4)];

insetAx3 = axes('Parent', fig3, 'Position', insetAbs, 'Color', 'w');
plt_single_MFD(37);

% === Second inset ===
insetRel = [0.5, 0.15, 0.25, 0.3];
insetAbs = [ ...
    mainPos(1) + insetRel(1) * mainPos(3), ...
    mainPos(2) + insetRel(2) * mainPos(4), ...
    insetRel(3) * mainPos(3), ...
    insetRel(4) * mainPos(4)];

insetAx4 = axes('Parent', fig3, 'Position', insetAbs, 'Color', 'w');
plt_single_MFD(8);



% === Helper functions ===


function [pl] = plt_MFD(M, color, label)
    dm = 0.1;
    m = 0 : dm : max(M);
    if max(m) < max(M)
        m = [m, max(M)];
    end
    for i=1:length(m)
        cn0(i)=length(find(M>=m(i)));
    end
   
    semilogy(m, cn0/cn0(1), 'color', color, 'linewidth', 2);
    hold on;
    pl = plot(m(end), cn0(end)/cn0(1), 'o', ...
    'Color', color, ...          % edge color
    'MarkerFaceColor', color, ...% fill color (same as edge)
    'MarkerSize', 12, ...        % larger marker
    'LineWidth', 2.5);           % line width (if connected lines exist)

    %scatter(m(end), cn0(end)/cn0(1), 150, color, 'filled', ...
     %   'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    text(m(end), cn0(end)/cn0(1), label, 'Color', 'w', 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
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

function [] = plt_mcorner_his()
    % === Load files ===
    dataDir = "./Data";
    file_name = {dir(fullfile(dataDir,"*.txt")).name};
    sites = cellfun(@(f) fileparts(f), file_name, "uni", 0);

    S = load("Case_Summary.mat");               % expects S.Case_Summary
    Case_Summary = S.Case_Summary;
    Data_Flag = load(fullfile(dataDir,"_Data_Flag.md"));  % assumes pure numeric ASCII

    % === Sort by col 12 (desc) ===
    [~, idx] = sort(Case_Summary(:,12), "descend");
    Case_Summary = Case_Summary(idx,:); 
    Data_Flag = Data_Flag(idx,:);
    if numel(sites) == numel(idx), sites = sites(idx); end

    % === Filter sci_flag ~= 1 (col 9 of Data_Flag) ===
    keep = Data_Flag(:,9) ~= 1;
    Case_Summary = Case_Summary(keep,:); 
    sites = sites(keep);

    % === Prep & plot ===
    bin_edges = -1:1:4; 
    bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2; %#ok<NASGU>
    P_LR = Case_Summary(:,12); 
    m_corner = Case_Summary(:,9);
    TGR_5 = P_LR <= 0.05; % GR_5 = P_LR > 0.05;

    % Plot histogram
    h = histogram(m_corner(TGR_5), bin_edges, ...
    'FaceColor', [0.85, 0.1, 0.1], 'EdgeColor', 'none', ...
    'FaceAlpha', 1);  % fully opaque

    ylim([0 8]); 
    xlim([-1 4]);
    xlabel('Estimated {\itm}_{corner}', 'Interpreter', 'tex'); 
    ylabel('Count');
    set(gca, 'FontSize', 14, 'Box', 'off', 'Color', 'none');

end




function []=plt_single_MFD(plt_idx)
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

    [~, sort_idx] = sort(Case_Summary(:,12), 'descend');
    sites = sites(sort_idx);
    Case_Summary = Case_Summary(sort_idx,:);


    site = sites(plt_idx);
    Case_Summary = Case_Summary(plt_idx,:);

        file = ['./Data/', site{1}, '.txt'];
        D = load(file); 
        D = sortrows(D, 1);
        D = D(D(:,1) > -1e4, :);
        T0 = D(:,1);
        M0 = D(:,5);
        i=1;
        mc = Case_Summary(i,1);
        b_GR = Case_Summary(i,3);
        b_TGR = Case_Summary(i,6);
        mcorner_TGR = Case_Summary(i,9)-mc;
        P_LR = Case_Summary(i,12);
        alpha_GR = Case_Summary(i,13);
        alpha_TGR = Case_Summary(i,14);
        Obs_max = Case_Summary(i,15);
        mode_max_GR = Case_Summary(i,16);
        max_GR_5 = Case_Summary(i,17);
        max_GR_95 = Case_Summary(i,18);
        mode_max_TGR = Case_Summary(i,19);
        max_TGR_5 = Case_Summary(i,20);
        max_TGR_95 = Case_Summary(i,21);
        color = custom_cmap(1 + floor(round(P_LR, 2) * 100), :);

        idx = M0 - mc >= -1e-11;
        Dm = M0(idx) - mc;

        plt_single_MFD_for_fig1_inset(M0, Dm, mc, b_GR, b_TGR, mcorner_TGR, P_LR, alpha_GR, alpha_TGR);
        grid off;
        title(site);

end