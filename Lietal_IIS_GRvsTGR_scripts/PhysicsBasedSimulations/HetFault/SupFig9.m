clc,clear
close all
% File: make_plots.m
% MATLAB rewrite of your Julia script

% ---------- Load data ----------
S = load('InputParameters.mat');
Fault_a            = S.Fault_a;
Fault_b            = S.Fault_b;
Fault_NormalStress = S.Fault_NormalStress;
FaultCenter        = S.FaultCenter;   % expected size: [N x 2]


FontSize = 16;

FaultCount = numel(Fault_a);
XYLength   = round(sqrt(FaultCount));

MinX = min(FaultCenter(:,1));
MaxX = max(FaultCenter(:,1));
MinY = min(FaultCenter(:,2));
MaxY = max(FaultCenter(:,2));

% Grid spacing (assumes uniform grid of XYLength x XYLength)
FaultGapX = (MaxX - MinX) / (XYLength - 1);
FaultGapY = (MaxY - MinY) / (XYLength - 1);  % use Y extent for Y spacing

% Preallocate grids
ContourX   = zeros(XYLength, XYLength);
ContourY   = zeros(XYLength, XYLength);
ContourPlot = zeros(XYLength, XYLength);

% Integer grid indices for each fault center
Xposition = floor((FaultCenter(:,1) - MinX) / FaultGapX) + 1;
Yposition = floor((FaultCenter(:,2) - MinY) / FaultGapY) + 1;

% Clamp just in case of numeric edge cases
Xposition = max(min(Xposition, XYLength), 1);
Yposition = max(min(Yposition, XYLength), 1);

% Build coordinate grids
xv = MinX + (0:XYLength-1) * FaultGapX;
yv = MinY + (0:XYLength-1) * FaultGapY;
[CX, CY] = meshgrid(xv, yv);
ContourX = CX;
ContourY = CY;

%% -------------- Plot 1: a - b ----------------
FaultValue = Fault_a - Fault_b;


% Place values into the grid
ContourPlot(:) = 0;
for k = 1:FaultCount
    ContourPlot(Yposition(k), Xposition(k)) = FaultValue(k); % note row=y, col=x
end

levels = -0.004:0.0001:0.001;

figure1=figure('Units','inches','Position',[1 1 5 5], ...
              'Color','w','Renderer','painters');
contourf(ContourX, ContourY, ContourPlot, levels, 'LineStyle','none');
% Use viridis if available, otherwise fall back
try
    colormap(viridis); %#ok<*UNRCH>
catch
    try, colormap('viridis'); catch, colormap('copper'); end
end
cbar = colorbar;
set(cbar, 'Ticks', -0.004:0.001:0.001);
xticks([-1000 -500 0 500 1000]);
yticks([-1000 -500 0 500 1000]);
xlabel('X (m)', 'FontSize', FontSize);
ylabel('Y (m)', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
ylabel(cbar, '{\ita}_{RS} - {\itb}_{RS}', 'FontSize', FontSize);
set(cbar, 'FontSize', FontSize);
set(gcf, 'Color', 'w');
axis equal;
%% -------------- Plot 2: Normal Stress (MPa) ----------------
FaultValue = Fault_NormalStress;

ContourPlot(:) = 0;
for k = 1:FaultCount
    ContourPlot(Yposition(k), Xposition(k)) = FaultValue(k);
end

figure2=figure('Units','inches','Position',[1 1 5 5], ...
              'Color','w','Renderer','painters');

contourf(ContourX, ContourY, ContourPlot/1e6, 50, 'LineStyle','none'); % Pa -> MPa
try
    colormap(viridis);
catch
    try, colormap('viridis'); catch, colormap('copper'); end
end
cbar = colorbar;
set(cbar, 'Ticks', 0:4:20);
xticks([-1000 -500 0 500 1000]);
yticks([-1000 -500 0 500 1000]);
xlabel('X (m)', 'FontSize', FontSize);
ylabel('Y (m)', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
ylabel(cbar, 'Initial Normal Stress \sigma_0 (MPa)', 'FontSize', FontSize);
set(cbar, 'FontSize', FontSize);
set(gcf, 'Color', 'w');
axis equal;
%%
figure3=figure('Units','inches','Position',[1 1 5 5], ...
              'Color','w','Renderer','painters');
FaultValue = Fault_NormalStress.*(Fault_a - Fault_b);

ContourPlot(:) = 0;
for k = 1:FaultCount
    ContourPlot(Yposition(k), Xposition(k)) = FaultValue(k);
end

contourf(ContourX, ContourY, ContourPlot, 50, 'LineStyle','none'); % Pa -> MPa
try
    colormap(viridis);
catch
    try, colormap('viridis'); catch, colormap('copper'); end
end
cbar = colorbar;
set(cbar, 'Ticks', -8e4:2e4:0);
xticks([-1000 -500 0 500 1000]);
yticks([-1000 -500 0 500 1000]);
xlabel('X (m)', 'FontSize', FontSize);
ylabel('Y (m)', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
ylabel(cbar, '({\ita}_{RS} - {\itb}_{RS})\sigma_0 (Pa)', 'FontSize', FontSize);
set(cbar, 'FontSize', FontSize);
set(gcf, 'Color', 'w');
axis equal;
