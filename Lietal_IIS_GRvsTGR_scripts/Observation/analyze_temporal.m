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

% === Load summary and data flags ===
load('Case_Summary.mat');
Data_Flag = load('./Data/_Data_Flag.md');  % Ensure this file is formatted as a numeric matrix

% === Sort by column 12 of Case_Summary ===
[~, sort_idx] = sort(Case_Summary(:,12), 'descend');
sites = sites(sort_idx);
Case_Summary = Case_Summary(sort_idx, :);
Data_Flag = Data_Flag(sort_idx, :);

% === Filter by catalog size threshold ===
catalogsize = Case_Summary(:,2);
jkf = find(catalogsize > 6e2);  % Keep only those with catalog size > 600
sites = sites(jkf);
Case_Summary = Case_Summary(jkf, :);
Data_Flag = Data_Flag(jkf, :);
Sample_number = length(jkf);

% === Start parallel pool ===
parpool(20)

% === Analysis loop ===
parfor i = 1:Sample_number
    site_ind(i) = jkf(i);
    mc(i) = Case_Summary(i, 1);
    global_bGR(i) = Case_Summary(i, 3);
    global_bTGR(i) = Case_Summary(i, 6);
    global_mcorner_TGR(i) = Case_Summary(i, 9);
    time_flag(i) = Data_Flag(i, 6);
    sci_flag(i) = Data_Flag(i, 9);

    site = sites{i};
    file = ['./Data/', site, '.txt'];
    D = load(file);

    % === Run temporal analysis ===
    [t_complete, m_complete, median_PLR, median_bGR, median_bTGR, median_mcorner, ...
        median_EQrate_aGR, median_COV, median_LV, mmax_pre, mmax_obs] = ...
        case_analysis_temporal(D, mc(i), global_bGR(i), global_bTGR(i), global_mcorner_TGR(i), time_flag(i));

    % === Save figure ===
    sgtitle(site);
    savefig(['./Figure/', site, '.fig']);

    % === Normalize event order ===
    order_complete = (1:length(t_complete));
    order_complete = (order_complete - min(order_complete)) / (max(order_complete) - min(order_complete));
    [~, ind] = max(m_complete);
    order_complete = order_complete - order_complete(ind);

    % === Store outputs ===
    allt_complete{i} = t_complete;
    allorder_complete{i} = order_complete;
    allm_complete{i} = m_complete;

    allMedian_PLR{i} = median_PLR;
    allMedian_bGR{i} = median_bGR;
    allMedian_bTGR{i} = median_bTGR;
    allMedian_mcorner{i} = median_mcorner;
    allmedian_EQrate_aGR{i} = median_EQrate_aGR;
    allmedian_COV{i} = median_COV;
    allmedian_LV{i} = median_LV;

    allmedian_mmax_obs{i} = mmax_obs;
    allmedian_mmax_GR_mode_st1{i} = mmax_pre.GR_mode_st1;
    allmedian_mmax_GR_5_st1{i} = mmax_pre.GR_5_st1;
    allmedian_mmax_GR_95_st1{i} = mmax_pre.GR_95_st1;
    allmedian_mmax_GR_mode_st2{i} = mmax_pre.GR_mode_st2;
    allmedian_mmax_GR_5_st2{i} = mmax_pre.GR_5_st2;
    allmedian_mmax_GR_95_st2{i} = mmax_pre.GR_95_st2;
    allmedian_mmax_TGR_mode_st1{i} = mmax_pre.TGR_mode_st1;
    allmedian_mmax_TGR_5_st1{i} = mmax_pre.TGR_5_st1;
    allmedian_mmax_TGR_95_st1{i} = mmax_pre.TGR_95_st1;
    allmedian_mmax_TGR_mode_st2{i} = mmax_pre.TGR_mode_st2;
    allmedian_mmax_TGR_5_st2{i} = mmax_pre.TGR_5_st2;
    allmedian_mmax_TGR_95_st2{i} = mmax_pre.TGR_95_st2;

    % === Correlation calculation ===
    variable_matrix = [
        median_PLR', median_bGR', median_EQrate_aGR', ...
        median_bTGR', median_mcorner', median_COV', median_LV'];

    corr_row = zeros(1, size(variable_matrix, 2));  % Preallocate local row
    for k = 1:size(variable_matrix, 2)
        corr_row(k) = corr(variable_matrix(:,1), variable_matrix(:,k));
    end
    corrr(i, :) = corr_row;  % Assign once outside the inner loop

end

% === Shutdown parallel pool ===
delete(gcp('nocreate'))

% === Save outputs ===
save('OUTPUT_temporal.mat');
