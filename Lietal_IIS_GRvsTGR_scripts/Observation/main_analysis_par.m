clc;
clear;
close all;

namelist = dir('./Data/*.txt');
file_name = {namelist.name};
case_num = length(file_name);
sites = cell(case_num, 1);
for i = 1:case_num
    [~, sites{i}] = fileparts(file_name{i});
end
Data_Flag=load('./Data/_Data_Flag.md');

lowest_mcs =Data_Flag(:,8);
Case_Summary=[];
parpool(40);
for i = 1:case_num
    i
    site = sites{i};  % sliced variable, allowed
    file = ['./Data/', site, '.txt'];
    D = load(file); 
    D = sortrows(D,1);
    jkf = find(D(:,1) > -1e4);
    D = D(jkf,:);
    T0 = D(:,1);
    M0 = D(:,5);
    Case_Summary(i,:) = analyze_each_case(T0, M0, lowest_mcs(i),site);
end
delete(gcp('nocreate'));
save('Case_Summary.mat', 'Case_Summary')
%% Local functions

function [OUTPUT]=analyze_each_case(T0,M0,lowest_mc,site)
%% Create figure
    figure('units','normalized','position',[0.05,0.05,1.8,0.9]);
    sgtitle(site);

    % Subplot 1: Original magnitude vs time
    subplot(1,4,1)
    scatter(T0, M0, 'filled');
    xlabel('Time');
    ylabel('M');
    grid on; box on; grid minor;
    set(gca, 'fontsize', 16);

    %% Magnitude processing
    diffm = diff(M0);
    nonzero_elements = diffm(diffm ~= 0);
    delta = min(abs(nonzero_elements)) / 2;

    dm = 0.1;
    m = floor(min(M0)/dm)*dm : dm : max(M0);
    n0 = hist(M0, m);
    [~, ind] = max(n0);
    Mc_MLE = m(ind) + 0.2;
    [Mc_MBS, b_es, number] = MBS_MLE_discrete(M0, delta);
    mc = max([Mc_MLE, Mc_MBS, lowest_mc]);
    

    idx = M0 - mc >= -1e-11;
    Dm = M0(idx) - mc;
    Dt = T0(idx);
    Num_above_mc=length(Dt);

    % Subplot 2: Processed magnitudes
    subplot(1,4,2)
    scatter(1:length(M0(idx)), M0(idx), 'filled');
    xlabel('Index');
    ylabel('M');
    grid on; box on; grid minor;
    set(gca, 'fontsize', 16);

    %% GR and TGR estimation
    L = length(Dm);
    mmin = 0;
    [b_GR, loglikelihood_GR] = Estimation_GR_discrete(Dm, mmin, delta);
    [b_TGR, mcorner_TGR, loglikelihood_TGR] = Estimation_TGR_gridsearch_discrete(Dm, mmin, delta, b_GR);
    R = loglikelihood_TGR - loglikelihood_GR;
    P_LR = 1 - chi2cdf(2 * R, 1);
    
    boostrap=1e3;
    max_GR_all=[];
    max_TGR_all=[];
    for k=1:boostrap
        Dm_rand = Dm(randi(L,1,L));
        [b_GR_bt(k),loglikelihood_GR]=Estimation_GR_discrete(Dm_rand,mmin,delta);
        [b_TGR_bt(k),mcorner_TGR_bt(k),loglikelihood_TGR]=Estimation_TGR_gridsearch_discrete(Dm_rand,mmin,delta,b_GR_bt(k));
        [max_GR,max_TGR]=get_mmax(L,b_GR_bt(k),b_TGR_bt(k),mcorner_TGR_bt(k),delta);
        max_GR_all=[max_GR_all,max_GR];
        max_TGR_all=[max_TGR_all,max_TGR];        
    end
    b_GR_percentiles = prctile(b_GR_bt, [5, 95]);
    b_TGR_percentiles = prctile(b_TGR_bt, [5, 95]);
    mcorner_TGR_percentiles = prctile(mcorner_TGR_bt, [5, 95]);
    OUTPUT_Parameter=[mc, Num_above_mc,b_GR,  b_GR_percentiles,b_TGR,b_TGR_percentiles, mcorner_TGR+mc, mc+mcorner_TGR_percentiles];
   
    Obs_max = max(Dm);
    mode_max_GR = log10(L) / b_GR;
    mode_max_TGR = get_mode_mmax_tgr(b_TGR, mcorner_TGR, L);
    mmax_GR_percentiles=prctile(max_GR_all, [5, 95]);
    mmax_TGR_percentiles=prctile(max_TGR_all, [5, 95]);
    
%     [max_GR_5, max_GR_95, max_TGR_5, max_TGR_95]=get_mmax_distribution(max_GR_all,max_TGR_all);
%     OUTPUT_mmax = mc + [Obs_max,mode_max_GR, max_GR_5, max_GR_95, mode_max_TGR, max_TGR_5, max_TGR_95];

OUTPUT_mmax = mc + [Obs_max,mode_max_GR, mmax_GR_percentiles, mode_max_TGR, mmax_TGR_percentiles];

    

    
    %% Bootstrap testing
    [alpha_GR,alpha_TGR] = get_alpha(L,b_GR,b_TGR,mcorner_TGR,P_LR,delta);
    OUTPUT_Test=[P_LR, alpha_GR,alpha_TGR];


    % Subplot 3: MFD
    subplot(1,4,3)
    plt_single_MFD_for_fig1(M0, Dm, mc, b_GR, b_TGR, mcorner_TGR, P_LR, alpha_GR, alpha_TGR);

    % Subplot 4: Predicted vs Observed mmax
    subplot(1,4,4)
    maglim_min=-4;
    maglim_max=6;
    plot([maglim_min maglim_max], [maglim_min maglim_max], '--', 'linewidth', 5, 'color', [0.75 0.75 0.75]); hold on;

    h(1) = plot([OUTPUT_mmax(3), OUTPUT_mmax(4)], [OUTPUT_mmax(1), OUTPUT_mmax(1)], ...
                'color', 1/255*[0 139 200], 'linewidth', 2.5);
    h(2) = scatter(OUTPUT_mmax(2), OUTPUT_mmax(1), 50, 1/255*[0 139 200], 'filled');
    h(3) = plot([OUTPUT_mmax(6), OUTPUT_mmax(7)], [OUTPUT_mmax(1), OUTPUT_mmax(1)], ...
                'color', 1/255*[178 91 0], 'linewidth', 1.5);
    h(4) = scatter(OUTPUT_mmax(5), OUTPUT_mmax(1), 30, 1/255*[178 91 0], 'filled');

    axis equal;
    xlim([maglim_min maglim_max]);
    ylim([maglim_min maglim_max]);
    xticks(maglim_min:maglim_max);
    yticks(maglim_min:maglim_max);
    %legend(h, 'GR 90%', 'GR Mode', 'TGR 90%', 'TGR Mode');
    grid on; box on; grid minor;
    xlabel('Predicted {\it m}_{max}');
    ylabel('Observed {\it m}_{max}');
    set(gca, 'fontsize', 16);
    saveas(gcf, fullfile('./Figure', [site '.png']));
    OUTPUT=[OUTPUT_Parameter,OUTPUT_Test,OUTPUT_mmax];
end
function Dm = gentgr(beta, L, Mmin, M_corner)
    U10 = rand(1, L);
    U20 = rand(1, L);
    U1 = -M_corner * log(U10);
    U2 = Mmin * U20.^(-1/beta);
    DM = min(U1, U2);
    Dm = 2/3 * log10(DM) -  6.07;
end

function [mode_mmax] = get_mode_mmax_tgr(b_TGR, mcorner_TGR, L)
    M_corner = 10^(1.5 * mcorner_TGR + 9.1);
    mmin_set = 0;
    mmin_true = mmin_set;
    Mmin = 10^(1.5 * mmin_true + 9.1);

    beta = 2/3 * b_TGR;
    mtest = 0:0.01:log10(L)/b_TGR;
    M = 10.^(1.5 * mtest + 9.1);
    misfit = abs((Mmin ./ M).^beta .* exp((Mmin - M) ./ M_corner) - 1/L);
    [~, ind] = min(misfit);
    mode_mmax = mtest(ind);
end

function [max_GR,max_TGR]=get_mmax(L,b_GR,b_TGR,mcorner_TGR,delta)
    rand_num = 1e3;
    max_GR = zeros(1, rand_num);
    max_TGR = zeros(1, rand_num);
    parfor k = 1:rand_num
        mmin_set = 0;
        mmin_true = mmin_set - delta;
        Mmin_set = 10^(1.5 * mmin_true + 9.1);

        % GR synthetic catalog
        Dm_GR = -1 / b_GR * log10(rand(1, L)) + mmin_true;
        Dm_GR = roundn(Dm_GR, round(log10(2 * delta)));
        max_GR(k) = max(Dm_GR);

        % TGR synthetic catalog
        M_corner = 10^(1.5 * mcorner_TGR + 9.1);
        Dm_TGR = gentgr(2/3 * b_TGR, L, Mmin_set, M_corner);
        Dm_TGR = roundn(Dm_TGR, round(log10(2 * delta)));
        max_TGR(k) = max(Dm_TGR);
    end
end

    
function [alpha_GR,alpha_TGR]=get_alpha(L,b_GR,b_TGR,mcorner_TGR,P_LR,delta)
    rand_num = 3e3;
    P_LR_rand_GR = zeros(1, rand_num);
    P_LR_rand_TGR = zeros(1, rand_num);
    parfor k = 1:rand_num
        mmin_set = 0;
        mmin_true = mmin_set - delta;
        Mmin_set = 10^(1.5 * mmin_true + 9.1);

        % GR synthetic catalog
        Dm_GR = -1 / b_GR * log10(rand(1, L)) + mmin_true;
        Dm_GR = roundn(Dm_GR, round(log10(2 * delta)));
        [b_GR_rand_est, loglikelihood_GR] = Estimation_GR_discrete(Dm_GR, mmin_set, delta);
        [b_TGR_rand_est, mcorner_TGR_rand_est, loglikelihood_TGR] = Estimation_TGR_gridsearch_discrete(Dm_GR, mmin_set, delta, b_GR_rand_est);
        R = loglikelihood_TGR - loglikelihood_GR;
        P_LR_rand_GR(k) = 1 - chi2cdf(2 * R, 1);

        % TGR synthetic catalog
        M_corner = 10^(1.5 * mcorner_TGR + 9.1);
        Dm_TGR = gentgr(2/3 * b_TGR, L, Mmin_set, M_corner);
        Dm_TGR = roundn(Dm_TGR, round(log10(2 * delta)));
        [b_GR_rand_est, loglikelihood_GR] = Estimation_GR_discrete(Dm_TGR, mmin_set, delta);
        [b_TGR_rand_est, mcorner_TGR_rand_est, loglikelihood_TGR] = Estimation_TGR_gridsearch_discrete(Dm_TGR, mmin_set, delta, b_GR_rand_est);
        R = loglikelihood_TGR - loglikelihood_GR;
        P_LR_rand_TGR(k) = 1 - chi2cdf(2 * R, 1);
    end

    alpha_GR = mean(P_LR_rand_GR < P_LR);
    alpha_TGR = mean(P_LR_rand_TGR > P_LR);
    
end