
clc;
clear;

% figure('units','normalized','position',[0.1,0.1,0.9,0.6])

%% Load file names
namelist1 = dir('*_out.txt');
file_name1 = {namelist1.name};
l = length(file_name1); 
z = 0;

namelist2 = dir('*_nihe.txt');
file_name2 = {namelist2.name};

Sample_number = l;

%% Create colormap
color = nan(Sample_number, 3);
color_end = 1/255 * [0, 0, 255];
color_begin = 1/255 * [135, 206, 235];

for i = 1:Sample_number
    color(i, :) = color_begin + (color_end - color_begin) * i / Sample_number;
end

rank_number = [4; 2; 8; 7; 9; 6; 3; 5; 1]; % Based on sample number

%% Plot data (log-log plot)
subplot(1, 2, 1);

for j = 1:l
    i = find(rank_number == j);
    D = load(file_name2{i});
    
    w1(i) = D(1);
    w2(i) = D(2);
    depth_min(i) = D(3);
    depth_max(i) = D(4);
    
    xi = linspace(depth_min(i), depth_max(i), 20);
    yi = w1(i) .* xi .^ (w2(i));
    
    h(j) = loglog(xi/1e6, yi, 'LineWidth', 2.5, 'Color', color(rank_number(i), :));
    hold on;
end

str = [
    "Leurer and Dvorkin, 2006";
    "Domenico, 1977";
    "Holt et al., 2005";
    "Linneman et al., 2021";
    "Rajabi and Sharifipour, 2017";
    "Chien and Oh, 2000";
    "Mao et al., 2022";
    "Debnath et al., 2022";
    "Zimmer et al., 2007"
];

legend(h, str);

xlabel('Confining Pressure (MPa)');
ylabel('\eta (Pa^{-1})');
xlim([5e-2, 1e3]);
ylim([1e-11, 1e-5]);
grid on;
box on;
set(gca, 'FontSize', 16);

%% Plot histograms
fig3 = figure('units', 'normalized', 'position', [0.1, 0.1, 0.6, 0.3]);

% First subplot: histogram of log10(a)
subplot(1, 2, 1);
histogram(log10(w1 .* 10.^(6 * w2)), 'FaceColor', 'b', 'BinWidth', 0.5);
xlabel('log_{10} \it{a}');
ylabel('Count');
xlim([-9, -6]);
ylim([0, 15]);
grid on;
box on;
grid minor;
set(gca, 'FontSize', 16);

% Second subplot: histogram of b
subplot(1, 2, 2);
histogram(w2, 'FaceColor', 'b', 'BinWidth', 0.5);
xlabel('\it{b}');
ylabel('Count');
xlim([-3, 0]);
ylim([0, 15]);
grid on;
box on;
grid minor;
set(gca, 'FontSize', 16);
