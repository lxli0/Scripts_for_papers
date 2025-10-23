clc; clear;
close all
load('OUTPUT_temporal.mat')
color = flip(turbo(Sample_number));  % Retain original color scheme

%% Categorization
cat_GR = [];
cat_TGR = [];



%% Categorization

for i = 1:Sample_number
    Temp_plr = allMedian_PLR{i};
    frac_below_005 = sum(Temp_plr < 0.05) / numel(Temp_plr);
    frac_below_010 = sum(Temp_plr < 0.10) / numel(Temp_plr);

    if frac_below_005 ==0
        cat_GR = [cat_GR, i];
    elseif frac_below_010==1
        cat_TGR = [cat_TGR, i];
    end
end
cat_TGR=[cat_TGR,17,14,10,16,13,15];
cat_tran=[3 4 2 12 1 11];

cat_GR=sort(cat_GR);
cat_TGR=sort(cat_TGR);
cat_tran=sort(cat_tran);


[~, TGR_order]  = sort(cellfun(@length, sites(cat_TGR)));
[~, tran_order] = sort(cellfun(@length, sites(cat_tran)));
[~, GR_order]   = sort(cellfun(@length, sites(cat_GR)));

cat_TGR_sorted  = cat_TGR(TGR_order);
cat_tran_sorted = cat_tran(tran_order);
cat_GR_sorted   = cat_GR(GR_order);

% Concatenate into final reordered index list
reorder_indices = [cat_TGR_sorted, cat_tran_sorted, cat_GR_sorted];

allMedian_PLR = allMedian_PLR(reorder_indices);
allMedian_bGR = allMedian_bGR(reorder_indices);
allMedian_bTGR = allMedian_bTGR(reorder_indices);
allMedian_mcorner = allMedian_mcorner(reorder_indices);

allt_complete = allt_complete(reorder_indices);
allorder_complete = allorder_complete(reorder_indices);
allm_complete = allm_complete(reorder_indices);

global_bGR = global_bGR(reorder_indices);
global_bTGR = global_bTGR(reorder_indices);
global_mcorner_TGR = global_mcorner_TGR(reorder_indices);
mc = mc(reorder_indices);

allmedian_mmax_obs=allmedian_mmax_obs(reorder_indices);
allmedian_mmax_GR_mode_st1=allmedian_mmax_GR_mode_st1(reorder_indices);
allmedian_mmax_GR_5_st1=allmedian_mmax_GR_5_st1(reorder_indices);
allmedian_mmax_GR_95_st1=allmedian_mmax_GR_95_st1(reorder_indices);
allmedian_mmax_GR_mode_st2=allmedian_mmax_GR_mode_st2(reorder_indices);
allmedian_mmax_GR_5_st2=allmedian_mmax_GR_5_st2(reorder_indices);
allmedian_mmax_GR_95_st2=allmedian_mmax_GR_95_st2(reorder_indices);
allmedian_mmax_TGR_mode_st1=allmedian_mmax_TGR_mode_st1(reorder_indices);
allmedian_mmax_TGR_5_st1=allmedian_mmax_TGR_5_st1(reorder_indices);
allmedian_mmax_TGR_95_st1=allmedian_mmax_TGR_95_st1(reorder_indices);
allmedian_mmax_TGR_mode_st2=allmedian_mmax_TGR_mode_st2(reorder_indices);
allmedian_mmax_TGR_5_st2=allmedian_mmax_TGR_5_st2(reorder_indices);
allmedian_mmax_TGR_95_st2=allmedian_mmax_TGR_95_st2(reorder_indices);


legend_labels = sites(reorder_indices);
num_each_cat = [length(cat_TGR), length(cat_tran), length(cat_GR)];
cum_num_each_cat = [0, length(cat_TGR), length(cat_TGR) + length(cat_tran)];

%% Figure 2 in maintext
figure('units','normalized','position',[0.1,0.1,0.9,0.65])

all_plot_handles = [];
all_plot_labels = {};
titles={'Persistent TGR';'GR–TGR Transition';'Persistent GR'};
for JJKK = 1:3
    plt_indices = cum_num_each_cat(JJKK) + (1:num_each_cat(JJKK));
    current_label = legend_labels(plt_indices);

    %% Panel positioning
    panel_width = 0.13;
    panel_height_top = 0.3;
    panel_height_bottom = 0.35;
    panel_bottom_top = 0.65;
    panel_bottom_bottom = 0.2;
    gap_x = 0.04;
    start_left = 0.08;

    % Left position for each category column
    left_pos = start_left + (JJKK - 1) * (panel_width * 2 + gap_x);

    %% === Row 1: p_LRT ===
    ax1 = axes('Position', [left_pos, panel_bottom_top, panel_width * 2, panel_height_top]); hold on;

    fill([-1, 1, 1, -1], [0, 0, 0.05, 0.05], [0.45 0.45 0.45], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    fill([-1, 1, 1, -1], [0.05, 0.05, 0.1, 0.1], [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

    for i = 1:num_each_cat(JJKK)
        idx = plt_indices(i);
        h = plot(allorder_complete{idx}, allMedian_PLR{idx}, '-', 'Color', color(idx,:), 'LineWidth', 2);
        all_plot_handles(end+1) = h;
        all_plot_labels{end+1} = current_label{i};
    end

    ylabel('{\it p}_{LRT}', 'FontName', 'Arial', 'FontSize', 16);
    title(titles{JJKK}, 'FontName', 'Arial', 'FontSize', 16);
    xlim([-1 1]); ylim([0 1]);
    set(gca, 'FontSize', 16, 'FontName', 'Arial', 'Box', 'off', 'TickDir', 'out');
    grid off;
    xlabel('Normalized Event Index');

    %% === Row 2 Left: mmax GR ===
    ax2 = axes('Position', [left_pos, panel_bottom_bottom, panel_width, panel_height_bottom]); hold on;
    plt_mmax_confidence_GR_b1(-4:1:6);

    for i = 1:num_each_cat(JJKK)
        idx = plt_indices(i);
        plot(allmedian_mmax_GR_mode_st1{idx}, allmedian_mmax_obs{idx}, '-', 'Color', color(idx,:), 'LineWidth', 1.5);
        scatter(allmedian_mmax_GR_mode_st1{idx}, allmedian_mmax_obs{idx}, 30, color(idx,:), 'filled');
    end

    xlabel('Predicted {\itm}_{\rm max}^{GR}', 'FontName', 'Arial', 'FontSize', 16);
    ylabel('Observed {\itm}_{\rm max}', 'FontName', 'Arial', 'FontSize', 16);
    axis equal; xlim([-4 6]); ylim([-4 6]);
    set(ax2, 'XTick', -4:2:4);

    set(gca, 'FontSize', 16, 'FontName', 'Arial', 'Box', 'off', 'TickDir', 'out');
    grid off;

    %% === Row 2 Right: mmax TGR ===
    ax3 = axes('Position', [left_pos + panel_width, panel_bottom_bottom, panel_width, panel_height_bottom]); hold on;

    plot([-4 6], [-4 6], '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);

    for i = 1:num_each_cat(JJKK)
        idx = plt_indices(i);
        plot(allmedian_mmax_TGR_mode_st1{idx}, allmedian_mmax_obs{idx}, '-', 'Color', color(idx,:), 'LineWidth', 1.5);
        scatter(allmedian_mmax_TGR_mode_st1{idx}, allmedian_mmax_obs{idx}, 30, color(idx,:), 'filled');
    end

    xlabel('Predicted {\itm}_{\rm max}^{TGR}', 'FontName', 'Arial', 'FontSize', 16);
    ylabel(ax3, '');
    set(ax3, 'YTickLabel', []);
    axis equal; xlim([-4 6]); ylim([-4 6]);
    set(ax3, 'XTick', -2:2:6);
    set(gca, 'FontSize', 16, 'FontName', 'Arial', 'Box', 'off', 'TickDir', 'out');
    grid off;
    
    % === Local legend for this category ===
    leg_ax = axes('Position', [left_pos, 0.11, panel_width * 2, 0.05], 'Visible', 'off');
    legend(leg_ax, all_plot_handles(cum_num_each_cat(JJKK) +1:cum_num_each_cat(JJKK) +num_each_cat(JJKK)), ...
           all_plot_labels(cum_num_each_cat(JJKK) +1:cum_num_each_cat(JJKK) +num_each_cat(JJKK)), ...
           'FontSize', 16, 'Location', 'southoutside', ...
           'Orientation', 'horizontal','NumColumns', 2,  'Box', 'off');
    
end


%% MFD paramters
figure('units','normalized','position',[0.1,0.1,0.9,0.8])

% Define manual subplot positions (3 rows × 3 columns)
lefts = [0.08, 0.38, 0.68];       % x-position for 3 columns
bottoms = 0.05+[0.66, 0.39, 0.13];     % y-position for 3 rows (reduced vertical spacing)
width = 0.25;                     % width of each subplot
height = 0.23;                    % height of each subplot (slightly reduced)

for JJKK = 1:3
    h = gobjects(1, num_each_cat(JJKK));
    
    % --- b_GR subplot ---
    subplot_pos = [lefts(JJKK), bottoms(1), width, height];
    ax1 = subplot('Position', subplot_pos);
    hold on;
    
    for i = 1:num_each_cat(JJKK)
        plt_rank = cum_num_each_cat(JJKK) + i;
        h(i) = plot(allorder_complete{plt_rank}, allMedian_bGR{plt_rank}, '-', 'color', color(plt_rank,:), 'linewidth', 2.5);
    end
    ylabel('{\it b}_{GR}');
    xlim([-1 1]); ylim([0.8 2.2]);
    xticks(-1:0.2:1);
    set(gca, 'fontsize', 16, 'XTickLabel', [], 'XMinorTick', 'on', ...
        'YMinorTick', 'on', 'XMinorGrid', 'on', 'TickDir', 'out', 'YMinorGrid', 'on');
    grid off;
     title(titles{JJKK});
     
    % --- b_TGR subplot ---
    subplot_pos = [lefts(JJKK), bottoms(2), width, height];
    ax2 = subplot('Position', subplot_pos);
    hold on;
    for i = 1:num_each_cat(JJKK)
        plt_rank = cum_num_each_cat(JJKK) + i;
        plot(allorder_complete{plt_rank}, allMedian_bTGR{plt_rank}, '-', 'color', color(plt_rank,:), 'linewidth', 2.5);
    end
    ylabel('{\it b}_{TGR}');
    xlim([-1 1]); ylim([0.5 2]);
    xticks(-1:0.2:1);
    set(gca, 'fontsize', 16, 'XTickLabel', [], 'XMinorTick', 'on', ...
        'YMinorTick', 'on', 'XMinorGrid', 'on', 'TickDir', 'out', 'YMinorGrid', 'on');
    grid off;
    
    % --- m_corner subplot ---
    subplot_pos = [lefts(JJKK), bottoms(3), width, height];
    ax3 = subplot('Position', subplot_pos);
    hold on;
    for i = 1:num_each_cat(JJKK)
        plt_rank = cum_num_each_cat(JJKK) + i;
        plot(allorder_complete{plt_rank}, allMedian_mcorner{plt_rank}, '-', 'color', color(plt_rank,:), 'linewidth', 2.5);
    end
    ylabel('{\it m}_{corner}');
    xlim([-1 1]); ylim([-4 5]);
    xticks(-1:0.2:1);
    xlabel('Normalized Event Index');
    set(gca, 'fontsize', 16, 'XMinorTick', 'on', ...
        'YMinorTick', 'on', 'XMinorGrid', 'on', 'TickDir', 'out', 'YMinorGrid', 'on');
    grid off;
    
    % === Local legend for this category ===
    leg_ax = axes('Position', [lefts(JJKK), 0.08, panel_width * 2, 0.05], 'Visible', 'off');
    legend(leg_ax, all_plot_handles(cum_num_each_cat(JJKK) +1:cum_num_each_cat(JJKK) +num_each_cat(JJKK)), ...
           all_plot_labels(cum_num_each_cat(JJKK) +1:cum_num_each_cat(JJKK) +num_each_cat(JJKK)), ...
           'FontSize', 16, 'Location', 'southoutside', ...
           'Orientation', 'horizontal','NumColumns', 2,  'Box', 'off');
end



function [hh, legendEntries] = plt_mmax_confidence_GR_b1(m)
    b = 1;
    
    % Compute the confidence bounds
    mm5 = -1 ./ b .* log10(-log(0.05));
    mm95 = -1 ./ b .* log10(-log(0.95));

    % Create the boundary curves
    lowerBound = m + mm5;
    upperBound = m + mm95;

    % Plot the 1:1 reference line
    plot(m, m, '--', 'LineWidth', 2.5, 'Color', [0.75 0.75 0.75]);
    hold on;

    % Fill the region between the bounds with gray
    fill([m, fliplr(m)], [upperBound, fliplr(lowerBound)], [0.7 0.7 0.7], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.2);

    % Prepare output
    hh = gobjects(1,1);
    hh(1) = plot(NaN, NaN, 'k'); % placeholder if you need legendEntries
    
    % Set limits
    xlim([-4 6]);
    ylim([-4 6]);
    
    legendEntries = {}; % Modify as needed
end
