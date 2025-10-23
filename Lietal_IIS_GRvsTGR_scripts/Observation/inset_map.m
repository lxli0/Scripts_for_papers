function []=inset_map_large(sites,color_save)

    % Settings
    inputFile = 'induced_seismicity_data.txt';

    % Read data
    data = readtable(inputFile, 'FileType', 'text', 'Delimiter', '\t');
    usedData = data(data.usecase == 1, :);  % Only usecase == 1

    % Assign colors using jet colormap
    operations = {'HF', 'EGS', 'WFD', 'UGS', 'Scientific'};
    nOps = length(operations);

    % Plot base map
    ax = worldmap([-60 70], [-140 160]);  % Same region as GMT
    setm(ax, 'FontColor', 'none');         % Remove lat/lon labels
    setm(ax, 'MLabelLocation', []);        % Hide meridians
    setm(ax, 'PLabelLocation', []);        % Hide parallels

    % Land = gray
    land = shaperead('landareas.shp', 'UseGeoCoords', true);
    geoshow(land, 'FaceColor', [0.9 0.9 0.9],'EdgeColor','k');  % Land areas in gray

    % Coastlines for detail
%     load coastlines
%     geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'none');

    % Plot points per operation
    for i = 1:nOps
        op = operations{i};
        opData = usedData(strcmp(usedData.operation, op), :);
        if ~isempty(opData)
            
            location=opData.location;
            lat = opData.lat;
            lon = opData.lon;

            for jj=1:length(location)
                    siteIdx0 = strcmp(location{jj}, sites); % Optional: get the index
                    [~,siteIdx]= max(siteIdx0);
                    color =color_save(siteIdx,:);
                    scatterm(lat(jj), lon(jj), 80, 'filled', 'MarkerFaceColor', color);
            end
        end
    end
    box off;
end

