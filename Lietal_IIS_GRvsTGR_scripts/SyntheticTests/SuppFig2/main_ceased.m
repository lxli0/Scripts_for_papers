clc; clear; close all

% ---------- Parameters ----------
rng(12);                 % Reproducibility
b          = 1;          % GR b-value
N          = 1e4;        % Catalog size
mmcut      = 2.5;        % Threshold magnitude to stop at first exceedance
Ltrials    = 10;         % Number of synthetic catalogs per panel

% ---------- Figure ----------
figure('units','normalized','position',[0.1,0.1,0.75,0.35])

% Panel 1: up to first M > mmcut
subplot(1,2,1)
runExperiment(b, N, mmcut, 0, Ltrials, ...
    'GR distribution, up to first {\itm} > 2.5');

% Panel 2: include 100 events after first M > mmcut
subplot(1,2,2)
runExperiment(b, N, mmcut, 100, Ltrials, ...
    'GR distribution, up to 100 events after first {\itm} > 2.5');

% ============================================================
% ===================== Helper Functions =====================
% ============================================================

function runExperiment(b, N, mmcut, Nafter, L, panelTitle)
    colos = turbo(L);

    % Preallocate
    P_LR_bt       = zeros(L,1);
    mcorner_TGR   = zeros(L,1);
    leg           = cell(L,1);
    plot([0, log10(N)/b],[N,1],'color',[0.75 0.75 0.75],'linewidth',2.5);
    hold on;
    for i = 1:L
        Dm_cal = gentgr_ceased(b, N, mmcut, Nafter);
        hh(i)=plotgr(Dm_cal, colos(i,:));

        % Your estimators (assumed to be on path)
        [~, loglikelihood_GR]              = Estimation_GR(Dm_cal, 0);
        [~, mcorner_TGR(i), loglikelihood_TGR] = Estimation_TGR(Dm_cal, 0);

        % Likelihood-ratio test (df=1)
        R           = loglikelihood_TGR - loglikelihood_GR;
        P_LR_bt(i)  = 1 - chi2cdf(2*R, 1);

%         leg{i} = sprintf('{\\itm}_{corner}^{est} = %.1f, {\\itp}_{LRT} = %.2f', ...
%                  mcorner_TGR(i), P_LR_bt(i));

    leg{i} = sprintf('{\\itp}_{LRT} = %.2f', P_LR_bt(i));
    end

    legend(hh,leg, 'Location','best');
    grid on; box on;
    xlabel('{\itm}');
    ylabel('{\itN}(≥{\it m})');
    set(gca,'fontsize',16);
    set(gca,'Yscale','log');

    % Axis limits chosen to match your original intent
    ylim([1 1e4]);
    xlim([0 6]);
    title(panelTitle);
end


function [pl]=plotgr(Dm,color)
    Catalogsize=length(Dm);
    dm=0.1;
    m=0:dm:max(Dm);
    n0=hist(Dm,m); 
    cn0=[];
    cn0(1)=Catalogsize;
    for i=2:length(m)
        cn0(i)=Catalogsize-sum(n0(1:i-1)); 
    end
    pl=semilogy(m,cn0,'color',color,'linewidth',2);
    hold on;
end

function Dm = gentgr_ceased(b, L, mmcut, Nafter)
    % Generate GR magnitudes: P(M >= m) ~ 10^{-b m}
    % Inverse CDF sampling with base-10:
    %   U ~ U(0,1),  M = -log10(U)/b
    U  = rand(1, L);
    Dm0 = -log10(U)/b;

    % Find the first event exceeding mmcut
    jk = find(Dm0 >= mmcut, 1, 'first');

    if isempty(jk)
        % If none exceed mmcut, just return the whole catalog (or choose your policy)
        Dm = Dm0;
        return
    end

    % Include Nafter events after first exceedance (clamped to catalog length)
    idxEnd = min(jk + Nafter, L);
    Dm = Dm0(1:idxEnd);
end

