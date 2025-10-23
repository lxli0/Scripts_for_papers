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
%{
% === Loop over cases ===
for i = 1:case_num
    if mod(i-1, 8) == 0
        figure('units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.5]);
        t = tiledlayout(2, 4, 'Padding', 'none', 'TileSpacing', 'compact');
    end

    site = sites{i};
    file = ['./Data/', site, '.txt'];
    D = load(file); 
    D = sortrows(D, 1);
    D = D(D(:,1) > -1e4, :);
    T0 = D(:,1);
    M0 = D(:,5);

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

    nexttile
    plt_single_MFD_for_fig1(M0, Dm, mc, b_GR, b_TGR, mcorner_TGR, P_LR, alpha_GR, alpha_TGR);
    title(site, 'Color', color);
    
end
%}
%%

plt_idx=[4 8 10 24 31 32 36 37];
case_num = length(plt_idx);
sites = sites(plt_idx);
Case_Summary = Case_Summary(plt_idx,:);
% === Loop over cases ===
figure('units', 'normalized', 'position', [0.1, 0.1, 0.4, 0.9]);
t = tiledlayout(4, 2 , 'Padding', 'none', 'TileSpacing', 'compact');
for i = 1:case_num


    site = sites{i};
    file = ['./Data/', site, '.txt'];
    D = load(file); 
    D = sortrows(D, 1);
    D = D(D(:,1) > -1e4, :);
    T0 = D(:,1);
    M0 = D(:,5);

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

    nexttile
    plt_single_MFD(M0, Dm, mc, b_GR, b_TGR, mcorner_TGR, P_LR, alpha_GR, alpha_TGR, 1);
    grid off;
    title(site);
    
end
%}

% === Colormap function ===
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

