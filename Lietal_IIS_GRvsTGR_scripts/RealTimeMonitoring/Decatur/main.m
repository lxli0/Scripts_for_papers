clc;
clear;
close all;

%% Load and sort data
D = load('09-Decatur (UGS).txt');
D = sortrows(D, 1);
T0 = D(:, 1);
M_all = D(:, 5);
event_order = (1:length(T0))';

%% Compute smallest non-zero magnitude difference
diff_m = diff(M_all);
delta = min(abs(diff_m(diff_m ~= 0))) / 2;

%% Apply magnitude cutoff
mc = -0.6;
valid_idx = find(M_all - mc >= -1e-11);
M_all = M_all(valid_idx);

%% Moving window parameters
window_increment = 100;
start_bin = 200;
num_windows = floor((length(M_all) - start_bin) / window_increment) + 1;

%% Preallocate storage
forward_pLRT = zeros(1, num_windows);
backward_pLRT_all = cell(1, num_windows);
backward_x_all = cell(1, num_windows);

%% Main moving window loop
for i = 1:num_windows

    if i == num_windows
        end_idx = length(M_all);
    else
        end_idx = start_bin + (i-1)*window_increment;
    end
    
    M_window = M_all(1:end_idx);
    Dm = M_window - mc;
    mmin = 0;
    
    % Fit GR and TGR models
    [b_GR, loglik_GR] = Estimation_GR_discrete(Dm, mmin, delta);
    [b_TGR, mcorner_TGR, loglik_TGR] = Estimation_TGR_gridsearch_discrete(Dm, mmin, delta, b_GR);
    
    % Compute LRT p-value
    LR_stat = 2 * (loglik_TGR - loglik_GR);
    forward_pLRT(i) = 1 - chi2cdf(LR_stat, 1);
    
    % Backward sub-windows within current window
    backward_pLRT = zeros(1, i);
    backward_x = zeros(1, i);
    
    for k = 1:i
        if i == num_windows
            len_back = length(M_all);
        else
            len_back = start_bin + (k-1)*window_increment;
        end
        
        backward_x(k) = end_idx - len_back + 1;
        M_subwindow = M_all(end_idx - len_back + 1:end_idx);
        Dm = M_subwindow - mc;
        
        [b_GR, loglik_GR] = Estimation_GR_discrete(Dm, mmin, delta);
        [b_TGR, mcorner_TGR, loglik_TGR] = Estimation_TGR_gridsearch_discrete(Dm, mmin, delta, b_GR);
        
        LR_stat = 2 * (loglik_TGR - loglik_GR);
        backward_pLRT(k) = 1 - chi2cdf(LR_stat, 1);
    end
    
    backward_pLRT_all{i} = backward_pLRT;
    backward_x_all{i} = backward_x;
end


%% Figure: Magnitudes and pLRT
figure('Units','normalized','Position',[0.05 0.05 0.3 0.45]);
sgtitle('Decatur');

%-----------------------------
% Subplot 1: pLRT and Magnitudes
%-----------------------------
subplot(2,1,1)
hold on;
xlabel('Event index');
yyaxis right;
ylabel('{\itp}_{LRT}');
set(gca,'YColor','k');
xlim([0 1.5e3]);
ylim([0 1]);

% Preplot empty scatter for background
scatter(1:2.2e3, NaN(1,2.2e3), 10, [0.75 0.75 0.75], 'filled');

for i = 1:num_windows
    if i == num_windows
        end_idx = length(M_all);
    else
        end_idx = start_bin + (i-1)*window_increment;
    end
    x_forward = start_bin + (0:i-1)*window_increment;
    y_forward = forward_pLRT(1:i);
    y_backward = backward_pLRT_all{i};
    x_backward = backward_x_all{i};
    
    x_backward_s(i)=x_backward(1);
    x_forward_s(i)=x_forward(end);
    y_backward_s(i)=y_backward(1);
    y_forward_s(i)=y_forward(end);
    if i == num_windows
        yyaxis right
        h1=plot(x_forward_s(2:end), y_backward_s(1:end-1), '-<', 'LineWidth', 2, 'Color', 1/255*[227 162 6] );hold on;
        h2=plot(x_forward_s, y_forward_s, '->', 'LineWidth', 2.5, 'Color',  1/255*[66 90 4] );


        yyaxis left
        scatter(1:length(M_all), M_all, 60, [0.75 0.75 0.75], 'filled');
        ylabel('Magnitude');
        ylim([-0.5 1.5]);
    end
    
end
legend([h2,h1], 'Expanding window', 'Sliding window','NumColumns',2, 'Location', 'north','box','on');

set(gca,'FontSize',14,'TickDir','out','Box','off');
set(gca,'YColor','k');

%% Estimate maximum magnitude predictions over time
subplot(2,1,2)
scatter(1:length(M_all), M_all, 60, [0.75 0.75 0.75], 'filled');
ylabel('Magnitude');
ylim([-0.5 1.5]);
set(gca,'YColor','k');

hold on;
xlabel('Event index');
xlim([0 1.5e3]);
set(gca,'YColor','k');

for i = 1:num_windows
    if i == num_windows
        end_idx = length(M_all);
    else
        end_idx = start_bin + (i-1)*window_increment;
    end
    M_window = M_all(1:end_idx);
    Dm = M_window - mc;
    
    [b_GR(i), loglik_GR] = Estimation_GR_discrete(Dm, 0, delta);
    [b_TGR(i), mcorner_TGR(i), loglik_TGR] = Estimation_TGR_gridsearch_discrete(Dm, 0, delta, b_GR(i));
    
    if i > 1
        N = length(Dm);
        mmax_pre.GR_mode(i) = mc + log10(N)/b_GR(i-1);
        mmax_pre.GR_5(i) = mmax_pre.GR_mode(i) - log10(-log(0.05))/b_GR(i-1);
        mmax_pre.GR_95(i) = mmax_pre.GR_mode(i) - log10(-log(0.95))/b_GR(i-1);
        
        [mmax_pre.TGR_mode(i), mmax_pre.TGR_5(i), mmax_pre.TGR_95(i)] = ...
            get_mmax_tgr(b_TGR(i-1), mcorner_TGR(i-1), N, mc);
    end
    mmax_obs(i) = max(M_window);
    end_idx_save(i) = end_idx;
end

% Shaded area for CI
x_fill = [end_idx_save(2:end), fliplr(end_idx_save(2:end))];
y_fill = [mmax_pre.TGR_5(2:end), fliplr(mmax_pre.TGR_95(2:end))];
fill(x_fill, y_fill, [0.7 0.7 1], 'FaceAlpha',0.4, 'EdgeColor','none');

% Mode and observed lines
h1=plot(end_idx_save(2:end), mmax_pre.TGR_mode(2:end), '-ob', 'LineWidth',2.5);
h2=plot(end_idx_save, mmax_obs, '-or', 'LineWidth',2.5);
ylim([-0.5 1.5]);
set(gca,'YColor','k','FontSize',14,'TickDir','out','Box','off');
legend([h1,h2],'{\itm}_{max}^{TGR}','{\itm}_{max}^{Observed}','numcolumns',2,'location','south','box','on');



%%
figure('Units','normalized','Position',[0.05 0.05 0.3 0.22]);
scatter(1:length(M_all), M_all, 60, [0.75 0.75 0.75], 'filled');
ylabel('Magnitude');
ylim([-0.5 4]);
set(gca,'YColor','k');

hold on;
xlabel('Event Index');
xlim([0 1.5e3]);
set(gca,'YColor','k');
% Shaded area for CI
x_fill = [end_idx_save(2:end), fliplr(end_idx_save(2:end))];
y_fill = [mmax_pre.GR_5(2:end), fliplr(mmax_pre.GR_95(2:end))];
fill(x_fill, y_fill, [0.7 0.7 1], 'FaceAlpha',0.4, 'EdgeColor','none');

% Mode and observed lines
h1=plot(end_idx_save(2:end), mmax_pre.GR_mode(2:end), '-ob', 'LineWidth',2.5);
h2=plot(end_idx_save, mmax_obs, '-or', 'LineWidth',2.5);
ylim([-0.5 4]);
set(gca,'YColor','k','FontSize',14,'TickDir','out','Box','off');
legend([h1,h2],'{\itm}_{max}^{GR}','{\itm}_{max}^{Observed}','numcolumns',2,'location','south','box','on');
title('Decatur')


%% Supporting Functions

function [mmax_mode, m5, m95] = get_mmax_tgr(b_TGR, mcorner, N, Mc)
    M_corner = 10^(1.5 * mcorner + 9.1);
    Mmin = 10^(9.1);
    beta = (2/3) * b_TGR;
    m_grid = 0:0.01:log10(N)/b_TGR;
    M = 10.^(1.5*m_grid + 9.1);
    misfit = abs((Mmin./M).^beta .* exp((Mmin - M)/M_corner) - 1/N);
    [~, idx] = min(misfit);
    mmax_mode = m_grid(idx) + Mc;
    
    max_sim = arrayfun(@(~) max(gentgr(beta, N, Mmin, M_corner)) + Mc, 1:1000);
    pctiles = prctile(max_sim, [5 95]);
    m5 = pctiles(1);
    m95 = pctiles(2);
end

function Dm = gentgr(beta, N, Mmin, M_corner)
    U1 = -M_corner * log(rand(1,N));
    U2 = Mmin * rand(1,N).^(-1/beta);
    M = min(U1,U2);
    Dm = 2/3*log10(M) - 6.07;
end
