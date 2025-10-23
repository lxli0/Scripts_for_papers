
clc; clear; close all;
figure('Units','inches','Position',[1 1 5 5], ...
              'Color','w','Renderer','painters');

visualfault(0)

%% === Nature-style defaults (global) ===
set(groot, ...
  'defaultFigureColor','w', ...
  'defaultAxesBox','off', ...
  'defaultAxesTickDir','out', ...
  'defaultAxesTickLength',[.02 .02], ...
  'defaultAxesLineWidth',0.75, ...
  'defaultAxesFontName','Helvetica', ...
  'defaultAxesFontSize',16, ...
  'defaultAxesXMinorTick','on', ...
  'defaultAxesYMinorTick','on', ...
  'defaultLineLineWidth',1.25, ...
  'defaultLineMarkerSize',6, ...
  'defaultStemLineWidth',1.25, ...
  'defaultLegendBox','off', ...
  'defaultColorbarBox','off');

%% === Input Parameters ===
sigma = 40e6;          % Shear stress (Pa)
Vdyn = 1e-2;           % Dynamic threshold slip rate (m/s)
duration_max = 1;      % Max event separation (s)

%% === Load Bulk Parameters Once (shared, read-only) ===
bulk_raw = readcell('Input_BulkFaultGeometry.txt', 'Delimiter', '\t');
ShearMod     = bulk_raw{2, 2};
PoissonRatio = bulk_raw{2, 3};
R_Density    = bulk_raw{2, 4};

geom_line  = cell2mat(bulk_raw(4, :));
St_L       = geom_line(4);
Dip_L      = geom_line(5);
MaxLeng    = geom_line(18);
Cell_Area  = MaxLeng^2;
nx         = round(St_L / MaxLeng);
ny         = round(Dip_L / MaxLeng);
x_coords   = linspace(0, St_L, nx) - St_L/2;
y_coords   = linspace(0, Dip_L, ny) - Dip_L/2;

% Loop through all test names
processedFiles = { ...
  'Processed_1_1_1_flowrate0.5_a_barrier5.0e-03_dfi-1.0e-01_thetai1.0e+00.mat'; ...
  'Processed_1_1_3_flowrate0.5_a_barrier5.0e-03_dfi5.0e-02_thetai1.0e+00.mat' ...
};
colos = 1/255 * [
    71, 213, 135; 
    235, 0, 153;
];


%%
fig1 = figure('Units','inches','Position',[1 1 18 4], ...
              'Color','w','Renderer','painters');

for jkjk = 1:2
    processedFile = processedFiles{jkjk};
    load(processedFile, 'flow_rate', 'catalogs', 'V_max_save', 't_mean', ...
        'V_plot', 'slip_plot', 'hydrau_diff', 'step_range_plt', ...
        'x_coords', 'y_coords', 'title_str', ...
        'MergeTimeCriteria', 'MergeDistanceCriteria');
    
    EQ_CAT = catalogs{3}; % [time, mag, weightedX, weightedY, firstX, firstY, farX, farY]
    T_EQ = EQ_CAT(:,1) / 86400; % convert to days
    M_EQ = EQ_CAT(:,2);

    % Plot
    figure(fig1);
    subplot(1,3,jkjk)
    hh(1)=plot([1e1 4e3],[4.7236 4.7236],'--k','LineWidth', 2.5);hold on;
    hh(jkjk+1)=scatter(T_EQ, M_EQ, 80*1.^M_EQ, colos(jkjk,:),'LineWidth', 2.5);    

    ylabel('Magnitude');
    xlabel('Time (day)');
    ylim([0 5]);
    set(gca, 'XScale', 'log');
    xlim([1e1 4e3]);
    grid off; box off; 
    natureify(gca);
end

subplot(1,3,1)
leg = {'Tectonic','Understressed'};
lgd = legend(leg, 'Location', 'southeast', 'FontSize', 16, ...
              'Box', 'off', 'NumColumns', 1);
subplot(1,3,2)
leg = {'Tectonic','Overstressed'};
lgd = legend(leg, 'Location', 'southeast', 'FontSize', 16, ...
              'Box', 'off', 'NumColumns', 1);


subplot(1,3,3)
linesty= '-';
hh2(1)=plot_simple_FMD(4.7236*ones(1,100), linesty, 'k', 5);

for jkjk=1:2
    processedFile = processedFiles{jkjk};
    load(processedFile, 'flow_rate', 'catalogs', 'V_max_save', 't_mean', ...
        'V_plot', 'slip_plot', 'hydrau_diff', 'step_range_plt', ...
        'x_coords', 'y_coords', 'title_str', ...
        'MergeTimeCriteria', 'MergeDistanceCriteria');

    EQ_CAT = catalogs{3};  % [time, mag, weightedX, weightedY, firstX, firstY, farX, farY]
    T_EQ = EQ_CAT(:, 1) / 86400; 
    M_EQ = EQ_CAT(:, 2);

    X_EQ = zeros(size(EQ_CAT, 1), 3);
    Y_EQ = zeros(size(EQ_CAT, 1), 3);
    X_EQ(:, 1) = EQ_CAT(:, 3);  X_EQ(:, 2) = EQ_CAT(:, 5);  X_EQ(:, 3) = EQ_CAT(:, 7);
    Y_EQ(:, 1) = EQ_CAT(:, 4);  Y_EQ(:, 2) = EQ_CAT(:, 6);  Y_EQ(:, 3) = EQ_CAT(:, 8);

    colo   = colos(jkjk,:);
    
    
    hh2(jkjk+1)=plot_simple_FMD(M_EQ, linesty, colo,2.5); 
    xlim([1 5]); ylim([1/50 1]);
    natureify(gca);
end
leg = {'Tectonic','Understressed','Overstressed'};
legend(hh2,leg,'Location','southwest','FontSize',16,'Box','off');


%% ================= Helper Functions =================

function [pl] = plot_simple_FMD(M, linesty, colo, linwid)
    dm = 0.1;
    m = 0 : dm : max(M);
    if max(m) < max(M)
        m = [m, max(M)];
    end
    for i=1:length(m)
        cn0(i)=length(find(M>=m(i)));
    end
   
    pl = semilogy(m, cn0/cn0(1), linesty, 'LineWidth', linwid, 'Color', colo); hold on;
    
    hold on;
    xlabel('{\itm}'); ylabel('{\itP}(≥{\it m})');
end

function natureify(ax)
% Minimal Nature-style polish for one axes (and its colorbar if present)
if nargin==0 || isempty(ax), ax = gca; end
set(ax, 'TickDir','out', 'TickLength',[.02 .02], 'LineWidth',1.5, ...
        'Box','off', 'Layer','top', 'FontName','Helvetica', 'FontSize',16);
grid(ax,'off');
ax.LabelFontSizeMultiplier = 1.0;
ax.TitleFontSizeMultiplier = 1.0;
ax.XAxis.Exponent = 0; 
ax.YAxis.Exponent = 0; 

% Make any sibling colorbar match ticks-out and 16 pt font
cbs = findall(ax.Parent,'Type','ColorBar');
for cb = reshape(cbs,1,[])
    cb.TickDirection = 'out';
    cb.LineWidth = 0.75;
    cb.FontName = 'Helvetica';
    cb.FontSize = 16;
end
end
