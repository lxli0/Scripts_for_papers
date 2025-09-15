function [hh,Vtotalmax]=plt_Otaniemi_operation(color_ope)
    % File paths
    pump_folder = './Data_St1/';
    pump_filename = fullfile(pump_folder, 'st1_pump_data_final.xlsx');
    bleed_filename = fullfile(pump_folder,'./Bleed_offs.xlsx');

    % Load pump data
    if exist(pump_filename, 'file') == 2
        pump_dat = readtable(pump_filename, 'VariableNamingRule', 'preserve');
    else
        error('File not found: %s', pump_filename);
    end

    % Load bleed-off data
    if exist(bleed_filename, 'file') == 2
        bleed_dat = readtable(bleed_filename, 'VariableNamingRule', 'preserve');
    else
        error('File not found: %s', bleed_filename);
    end

    % Define reference start time
    ref_time = datetime(2018, 6, 4);

    % Extract time, cumulative injection volume, and pressure
    time = pump_dat.('DATE_TIME_HEL'); % Timestamp
    time = days(time - ref_time); % Convert to days since June 4, 2018
    Total_Volume = pump_dat.('CUM_VOL_m3'); % Cumulative injection volume

    % Ensure cumulative volume is monotonically increasing
    jkf = find(diff(Total_Volume) < 0);
    Total_Volume(jkf + 1) = nan;
    Total_Volume(jkf + 2) = nan;

    % Extract pressure data
    pres_organic = pump_dat.('INJ_PRESSURE_bars'); % Pressure in bars
    pres_mpa = pres_organic * 0.1; % Convert to MPa

    % Extract bleed-off data
    bleed_time = bleed_dat.('END_BLEED');  % Time of bleed-off
    bleed_time = days(bleed_time - ref_time); % Convert to days since June 4, 2018
    bleed_vol = bleed_dat.('VOL_BLEED_M3'); % Volume of bleed-off

    % Subtract bleed-off volumes at corresponding times
    for i = 1:length(bleed_time)
        % Find closest time index in pump data
        [~, idx] = min(abs(time - bleed_time(i)));

        % Subtract the bleed-off volume from total injected volume
        Total_Volume(idx:end) = Total_Volume(idx:end) - bleed_vol(i);
    end

    yyaxis left;
    hh(1)=plot(time, pres_mpa,'-', 'color',color_ope(1,:), 'LineWidth', 2);
    ylabel('Well-head Pressure (MPa)');
    grid on;
    ylim([0 100]);
    hold on;
    ax = gca; % Get current axes
    ax.YColor = 'k';
    xlabel('Days Since 2018/6/4');
    set(gca,'fontsize',16)

    % Plot Cumulative Injection Volume vs. Time (After Bleed-off Correction)
    yyaxis right;
    hh(2)=plot(time, Total_Volume/max(Total_Volume), '-','color',color_ope(2,:), 'LineWidth', 2.5);hold on;
    grid on;
    ax = gca; % Get current axes
    ax.YColor = 'k';
    Vtotalmax=max(Total_Volume);
end