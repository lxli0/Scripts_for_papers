function [] = visualfault(barrier_flag)
    % ---------- Load data ----------
    S = load('InputParameters.mat');
    Fault_a     = S.Fault_a;
    Fault_b     = S.Fault_b;
    FaultCenter = S.FaultCenter;

    % ---------- Compute a - b and downsample ----------
    FaultValue = Fault_a - Fault_b;
    ds = 4;
    idx = 1:ds:numel(FaultValue);
    FaultValue = FaultValue(idx);
    FaultCenter = FaultCenter(idx,:);

    % ---------- Determine grid ----------
    % Estimate grid size assuming roughly square grid
    XYLength = round(sqrt(numel(FaultValue)));
    % Create regular grid from scattered fault centers
    [xq, yq] = meshgrid(...
        linspace(min(FaultCenter(:,1)), max(FaultCenter(:,1)), XYLength), ...
        linspace(min(FaultCenter(:,2)), max(FaultCenter(:,2)), XYLength));
    Fq = griddata(FaultCenter(:,1), FaultCenter(:,2), FaultValue, xq, yq, 'natural');

    % Offset to center like your scatter code
    xq = xq + 1050;
    yq = yq + 1050;
    zq = zeros(size(Fq));

    hold on;

    % ---------- Plot fault surface (pcolor style) ----------
    h = surf(xq, yq, zq, Fq, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp', ...
        'FaceLighting', 'none');

    % Custom gray→black colormap (log-scaled spacing)
    cmap_len = 256;
    % log10 spacing from 0.99 → 0.5 to compress light range and expand dark
     cmap_log = logspace(log10(0.99), log10(0.1), cmap_len)'; 
    cmap = [cmap_log cmap_log cmap_log];  % grayscale
    colormap(cmap);
    caxis([min(FaultValue(:)) 0]);

    % ---------- Draw fault boundary ----------
    x_rect = [0, 2100, 2100, 0, 0];
    y_rect = [0, 0, 2100, 2100, 0];
    plot3(x_rect, y_rect, zeros(size(x_rect)), 'k', 'LineWidth', 1.5);

    % ---------- Injection well ----------
    wellRadius = 30;
    wellHeight = 1000;
    wellCenter = [1050, 1050, 0];
    [wellX, wellY, wellZ] = cylinder(wellRadius, 100);
    wellZ = wellZ * wellHeight;
    wellX = wellX + wellCenter(1);
    wellY = wellY + wellCenter(2);
    wellZ = wellZ + wellCenter(3);
    surf(wellX, wellY, wellZ, 'FaceColor', 'k', 'EdgeColor', 'none');

    % ---------- Barrier or boundary ----------
    if barrier_flag == 1
        for xi = 0:100:2000
            for yi = 0:100:2000
                if mod(xi,400) < 100 || mod(yi,400) < 100
                    fill3([xi, xi+100, xi+100, xi], [yi, yi, yi+100, yi+100], ...
                          [0, 0, 0, 0], [0.75 0.75 0.75], 'EdgeColor','none');
                end
            end
        end
    else
        plot3([-100 -100],[0 2100],[0 0],'k','LineWidth',1.5);
        plot3([-50 -150],[0 0],[0 0],'k','LineWidth',1.5);
        plot3([-50 -150],[2100 2100],[0 0],'k','LineWidth',1.5);
        plot3([0 2100],[-100 -100],[0 0],'k','LineWidth',1.5);
        plot3([0 0],[-50 -150],[0 0],'k','LineWidth',1.5);
        plot3([2100 2100],[-50 -150],[0 0],'k','LineWidth',1.5);
    end

    % ---------- Colorful concentric circles ----------
    numCircles = 10;
    circleRadii = linspace(wellRadius, wellRadius*10, numCircles);
    colorm = parula(numCircles);
    for i = 1:numCircles
        theta = linspace(0, 2*pi, 100);
        xCircle = wellCenter(1) + circleRadii(i) * cos(theta);
        yCircle = wellCenter(2) + circleRadii(i) * sin(theta);
        zCircle = zeros(size(xCircle)); % z=0 plane
        plot3(xCircle, yCircle, zCircle, 'Color', colorm(numCircles+1-i,:), 'LineWidth', 1.5);
    end

    % ---------- Blue downward arrow ----------
    arrowStart = wellCenter + [100, 0, wellHeight/2+100];
    arrowEnd = arrowStart + [0, 0, -500];
    quiver3(arrowStart(1), arrowStart(2), arrowStart(3), ...
            arrowEnd(1)-arrowStart(1), arrowEnd(2)-arrowStart(2), arrowEnd(3)-arrowStart(3), ...
            'w', 'LineWidth', 2, 'MaxHeadSize', 20);

    % ---------- Horizontal direction arrows ----------
    quiver3(100, 1950, 0, 400, 0, 0, 'w', 'LineWidth', 1.5, 'MaxHeadSize', 1.5);
    quiver3(500, 2250, 0, -400, 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 1.5);

    % ---------- Appearance ----------
    daspect([1 1 1]);
    axis tight; axis off;
    set(gcf, 'Color', 'w');

    % ---------- 3D View ----------
    view(315, 30);
    camproj('perspective');
    lighting gouraud;
    material dull;
    camlight headlight;
end

