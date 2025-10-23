
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
% Max plotting time
plt_time_max = 4000;

% Loop through all test names
processedFiles = { ...
  'Catalog_Result_NoFlow' ...
     'Catalog_Result_Flow1' ...
  'Catalog_Result_Flow2' ...
  'Catalog_Result_UniformDrawdown1' ...
  'Catalog_Result_UniformDrawdown2' ...
};
starting_times=[0 4767 5295 4767 5295];
colos = 1/255 * [
    0 0 0
    71, 213, 135; 
    235, 0, 153;
    71, 213, 135; 
    235, 0, 153;
];


%% 
fig1 = figure('Units','inches','Position',[1 1 18 4], ...
              'Color','w','Renderer','painters');

for jkjk = 1:3
    EQ_CAT=load(processedFiles{jkjk});
    
    T_EQ = EQ_CAT(:,1) / 86400 + starting_times(jkjk); % convert to days
     M_EQ = EQ_CAT(:,2);

    % Plot
    figure(fig1);
    ax=subplot(1,3,1);
    hh(jkjk)=scatter(T_EQ, M_EQ, 40*1.^M_EQ, colos(jkjk,:),'LineWidth', 1.5);    
    hold on;

    ylabel('Magnitude');
    xlabel('Time (day)');
    ylim([0 4]);
    grid off; box off; 
    natureify(gca);

end
leg = {'Tectonic', 'Earl. Inter.','Late Inter.'};
lgd = legend(hh, leg, 'Location', 'southeast', 'FontSize', 16, ...
              'Box', 'on', 'NumColumns', 3);

% Set the background and edge colors
lgd.Color = 'w';          % white background
lgd.EdgeColor = 'w';      % black edge (default) — or pick another color

pos = get(ax, 'Position');  % Get current [x y width height]
pos(3) = pos(3) + 0.12;      % Move the axes upward
pos(2) = pos(2) + 0.05;      % Move the axes upward
set(ax, 'Position', pos);   % Apply the new position

natureify(gca);


ax=subplot(1,3,3)
hh=[];
linst={'-','-','-','--','--'};
linwidth=[5 2.5 2.5 2.5 2.5];
for jkjk = 1:5
    EQ_CAT=load(processedFiles{jkjk});
    
    M_EQ = EQ_CAT(:,2);
    colo   = colos(jkjk,:);
    
    hh2(jkjk)=plot_simple_FMD(M_EQ, jkjk, colo,linst{jkjk},linwidth(jkjk)); 
    xlim([1 4]); %ylim([1/ 1]);
    natureify(gca);
end
leg = {'Tectonic', ...
       'Earl. Inter. (near field)', ...
       'Late Inter. (near field)', ...
       'Earl. Inter. (far field)', ...
       'Late Inter. (far field)'};
legend(hh2, leg, 'Location','southwest', 'FontSize',16, 'Box','off');

fig1 = figure('Units','inches','Position',[1 1 18 4], ...
              'Color','w','Renderer','painters');

for jkjk = 1:1
    EQ_CAT=load(processedFiles{jkjk});
    
    T_EQ = EQ_CAT(:,1) / 86400 + starting_times(jkjk); % convert to days
     M_EQ = EQ_CAT(:,2);

    figure(fig1);
    ax=subplot(1,3,1);
    hh(jkjk)=scatter(T_EQ, M_EQ, 40*1.^M_EQ, colos(jkjk,:),'LineWidth', 1.5);    
    hold on;

    ylabel('Magnitude');
    xlabel('Time (day)');
    ylim([0 4]);
    grid off; box off; 
    natureify(gca);

end
plot([4767 4767],[0 4],'-','Color',colos(2,:),'linewidth',2.5);
plot([5295 5295],[0 4],'-','Color',colos(3,:),'linewidth',2.5);

xlim([4.2e3 5.8e3]);
pos = get(ax, 'Position');  % Get current [x y width height]
pos(2) = pos(2) + 0.06;      % Move the axes upward
pos(4) = pos(4) - 0.08;      % Move the axes upward
pos(3) = pos(3) - 0.1;      % Move the axes upward
set(ax, 'Position', pos);   % Apply the new position

%% ================= Helper Functions =================

function [pl] = plot_simple_FMD(M, tec_leg, colo,linst,linwidth)
    dm = 0.1;
    m = 0 : dm : max(M);
    if tec_leg==1
        m=m(1:end-1);
    end
    for i=1:length(m)
        cn0(i)=length(find(M>=m(i)));
    end
   
    pl = semilogy(m, cn0/cn0(1), linst, 'LineWidth', linwidth, 'Color', colo); hold on;
    
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
