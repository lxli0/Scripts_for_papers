clc;
clear;
close all;

%% Parameters
bs = 0.6:0.4:1.4;
mcorners = 0:0.5:4;
mc = 0;
Ls = [1e2, 1e3, 1e4, 1e5, 1e6];
colors = autumn(length(mcorners));

% Store parameters for each (b, mcorner)
p_store = zeros(length(bs), length(mcorners), 9); % 3 p's x 3 fits
figure('units','normalized','position',[0.1,0.1,0.6,0.9])
%% Main computation loop
for mm = 1:length(bs)
    b = bs(mm);
    
    
    fortitle=['{\itb} = ',num2str(b)];
    
    for jjkk = 1:length(mcorners)
        mcorner = mcorners(jjkk);
        color = colors(jjkk, :);

        % Get statistics and fitting parameters
        [mode_max_TGR, max_TGR_5, max_TGR_95, p_mode, p_5, p_95] = compute_fitting_parameters(mc, b, mcorner, Ls);
        
        % Store all 9 parameters (3 for each fit)
        p_store(jjkk, mm, :) = [p_mode, p_5, p_95];

        % Plot results
        logLs_fit = min(log10(Ls)):0.5:max(log10(Ls));
    
    subplot(length(bs),3,(mm-1)*3+1);
    scatter(Ls, mode_max_TGR, 30, color, 'filled'); hold on;
    plot(10.^logLs_fit, p_mode(1) ./ (1 + exp(-p_mode(2) * (logLs_fit - p_mode(3)))), '-', 'Color', color,'linewidth',2);
    plot(10.^logLs_fit,logLs_fit./b,'-k','linewidth',2);
    set(gca, 'XScale', 'log'); xlabel('{\itN}({\itm}\geq0)');
    ylabel('{\itm}_{max}^{mode}'); set(gca, 'FontSize', 16); xlim([1e2 1e8]);
    grid on; box on; grid minor;
    xticks([1e2 1e3 1e4 1e5 1e6]);
    xlim([min(Ls) max(Ls)])
    
    subplot(length(bs),3,(mm-1)*3+2);
    scatter(Ls, max_TGR_5, 30, color, 'filled'); hold on;
    plot(10.^logLs_fit, p_5(1) ./ (1 + exp(-p_5(2) * (logLs_fit - p_5(3)))), '-', 'Color', color,'linewidth',2);
    plot(10.^logLs_fit,logLs_fit./b-1/b*log10(-log(0.05)),'-k','linewidth',2);
    set(gca, 'XScale', 'log'); xlabel('{\itN}({\itm}\geq0)');
    ylabel('{\itm}_{max}^{5%}'); set(gca, 'FontSize', 16); xlim([1e2 1e8]);
    grid on; box on; grid minor;
    title(fortitle);
    xticks([1e2 1e3 1e4 1e5 1e6]);
    xlim([min(Ls) max(Ls)])
    
    subplot(length(bs),3,(mm-1)*3+3);
    scatter(Ls, max_TGR_95, 30, color, 'filled'); hold on;
    plot(10.^logLs_fit, p_95(1) ./ (1 + exp(-p_95(2) * (logLs_fit - p_95(3)))), '-', 'Color', color,'linewidth',2);
    plot(10.^logLs_fit,logLs_fit./b-1/b*log10(-log(0.95)),'-k','linewidth',2);
    set(gca, 'XScale', 'log'); xlabel('{\itN}({\itm}\geq0)');
    ylabel('{\itm}_{max}^{95%}'); set(gca, 'FontSize', 16); xlim([1e2 1e8]);
    grid on; box on; grid minor;
    xticks([1e2 1e3 1e4 1e5 1e6]);
    xlim([min(Ls) max(Ls)])
    
    
    end
end

for i=1:length(bs)*3
    subplot(length(bs),3,i);
    ylim([0 5]);
end
% Add colorbar with discrete colors
N = length(mcorners);
colormap(autumn(N));  % Create colormap with N discrete colors
cb = colorbar('Position', [0.92 0.11 0.02 0.82]);  % Positioned as in first version

% Calculate exact colorband edges and centers
edges = linspace(min(mcorners), max(mcorners), N+1);  % N+1 edges for N color bands
centers = edges(1:end-1) + diff(edges)/2;  % Centers of each color band

% Set up colorbar properties
caxis([min(mcorners), max(mcorners)]);  % Map data range to full colorbar
cb.Ticks = centers;  % Place ticks at centers of each color band
cb.TickLabels = arrayfun(@num2str, mcorners, 'UniformOutput', false);  % Your labels
cb.Label.String = '{\itm}_{corner}';
cb.Label.FontSize = 16;

%% Function Definitions
function [mode_max_TGR, max_TGR_5, max_TGR_95, p_mode, p_5, p_95] = compute_fitting_parameters(mc, b, mcorner, Ls)
    % Preallocate result arrays
    num_L = length(Ls);
    mode_max_TGR = zeros(1, num_L);
    max_TGR_5 = zeros(1, num_L);
    max_TGR_95 = zeros(1, num_L);

    %% Compute statistics for each catalog size
    for i = 1:num_L
        L = Ls(i);
        [mode_max_TGR(i), max_TGR_5(i), max_TGR_95(i)] = ...
            get_mmax(mc, b, mcorner, L);
    end

    %% Fit logistic curve
    logLs = log10(Ls)';
    logistic_model = @(p, x) p(1) ./ (1 + exp(-p(2) * (x - p(3))));
    p0 = [4, 1, mean(logLs)];
    opts = optimset('Display', 'off');

    p_mode = lsqcurvefit(logistic_model, p0, logLs, mode_max_TGR', [], [], opts);
    p_5 = lsqcurvefit(logistic_model, p0, logLs, max_TGR_5', [], [], opts);
    p_95 = lsqcurvefit(logistic_model, p0, logLs, max_TGR_95', [], [], opts);
end   

function [mode_TGR, TGR_5, TGR_95] = get_mmax(mc, b, mcorner, L)
    rand_num = 1e3;
    mcorner_TGR = mcorner - mc;
    
    max_TGR = zeros(1, rand_num);

    for k = 1:rand_num
        mmin_true = 0;
        Mmin_set = 10^(1.5 * mmin_true + 9.1);

        % TGR Catalog
        M_corner = 10^(1.5 * mcorner_TGR + 9.1);
        Dm_TGR = gentgr(2/3 * b, L, Mmin_set, M_corner);
        max_TGR(k) = max(Dm_TGR);
    end

    % Compute TGR statistics
    sorted_TGR = sort(max_TGR);
    TGR_5 = mc + sorted_TGR(floor(rand_num * 0.05));
    TGR_95 = mc + sorted_TGR(floor(rand_num * 0.95));
    mode_TGR = mc + get_mode_mmax_tgr(b, mcorner_TGR, L);
end

function Dm = gentgr(beta, L, Mmin, M_corner)
    U10 = rand(1, L);
    U20 = rand(1, L);
    U1 = -M_corner * log(U10);
    U2 = Mmin * U20.^(-1 / beta);
    DM = min(U1, U2);
    Dm = 2/3 * log10(DM) - 6.07;
end

function mode_mmax = get_mode_mmax_tgr(b_TGR, mcorner_TGR, L)
    M_corner = 10^(1.5 * mcorner_TGR + 9.1);
    mmin_true = 0;
    Mmin = 10^(1.5 * mmin_true + 9.1);
    beta = 2/3 * b_TGR;

    mtest = 0:0.01:log10(L) / b_TGR;
    M = 10.^(1.5 * mtest + 9.1);
    misfit = abs((Mmin ./ M).^beta .* exp((Mmin - M) ./ M_corner) - 1 / L);

    [~, ind] = min(misfit);
    mode_mmax = mtest(ind);
end