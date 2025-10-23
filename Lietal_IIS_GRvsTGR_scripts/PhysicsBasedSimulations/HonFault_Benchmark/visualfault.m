 function []=visualfault(barrier_flag)
    % Define the corners of the rectangle for the horizontal fault plane in x-y
    x = [0, 2100, 2100, 0, 0]; % x varies
    y = [0, 0, 2100, 2100, 0]; % y varies
    z = [0, 0, 0, 0, 0]; % z is constant for the x-y plane

    % Parameters for the injection well
    wellRadius = 30; % Well radius is small compared to the rectangle
    wellHeight = 1000; % Same as the rectangle size
    wellCenter = [1050, 1050, 0]; % Center of the well at the rectangle center (in x-y plane)


    % Plot the rectangle in 3D space with a lighter grey color and transparency
    fill3(x, y, z, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Lighter grey color with 50% transparency
    hold on; % Hold on to plot more objects

    % Plot the dark grey areas without borders
    if barrier_flag==1
        for xi = 0:100:2000
            for yi = 0:100:2000
                if mod(xi, 400) < 100 || mod(yi, 400) < 100
                    % Fill with dark grey for specified areas, no borders
                    fill3([xi, xi+100, xi+100, xi], [yi, yi, yi+100, yi+100], [0, 0, 0, 0], [0.75 0.75 0.75], 'EdgeColor', 'none');
                end
            end
        end
        plot3([-100 -100],[0 100],[0 0],'k','linewidth',1.5);
        plot3([-50 -150],[0 0],[0 0],'k','linewidth',1.5);
        plot3([-50 -150],[100 100],[0 0],'k','linewidth',1.5);

        plot3([-100 -100],[900 1200],[0 0],'k','linewidth',1.5);
        plot3([-50 -150],[900 900],[0 0],'k','linewidth',1.5);
        plot3([-50 -150],[1200 1200],[0 0],'k','linewidth',1.5);
        
        
        plot3([0 100],[-100 -100],[0 0],'k','linewidth',1.5);
        plot3([0 0],[-50 -150],[0 0],'k','linewidth',1.5);
        plot3([100 100],[-50 -150],[0 0],'k','linewidth',1.5);

        plot3([900 1200],[-100 -100],[0 0],'k','linewidth',1.5);
        plot3([900 900],[-50 -150],[0 0],'k','linewidth',1.5);
        plot3([1200 1200],[-50 -150],[0 0],'k','linewidth',1.5);
        

    else
        plot3([-100 -100],[0 2100],[0 0],'k','linewidth',1.5);
        plot3([-50 -150],[0 0],[0 0],'k','linewidth',1.5);
        plot3([-50 -150],[2100 2100],[0 0],'k','linewidth',1.5);
        plot3([0 2100],[-100 -100],[0 0],'k','linewidth',1.5);
        plot3([0 0],[-50 -150],[0 0],'k','linewidth',1.5);
        plot3([2100 2100],[-50 -150],[0 0],'k','linewidth',1.5);
    end
    %}

    % Draw the injection well as a cylinder
    [wellX, wellY, wellZ] = cylinder(wellRadius, 100); % 100 sides for a smooth cylinder
    wellZ = wellZ * wellHeight; % Scale the height of the cylinder
    wellX = wellX + wellCenter(1); % Shift the cylinder to the center of the rectangle
    wellY = wellY + wellCenter(2);
    wellZ = wellZ + wellCenter(3);
    surf(wellX, wellY, wellZ, 'FaceColor', 'k', 'EdgeColor', 'none'); % Black color

    % Plot concentric circles on the fault plane with gradient colors
    hold on;
    numCircles = 10;  % Number of concentric circles
    circleRadii = linspace(wellRadius, wellRadius * 10, numCircles);  % Radii of concentric circles
    colorm = parula(10);
    for i = 1:numCircles
        theta = linspace(0, 2*pi, 100);
        xCircle = wellCenter(1) + circleRadii(i) * cos(theta);
        yCircle = wellCenter(2) + circleRadii(i) * sin(theta);
        zCircle = wellCenter(3) * ones(size(xCircle)); % Keep z constant for the fault plane

        plot3(xCircle, yCircle, zCircle, 'Color', colorm(11-i,:), 'LineWidth', 1.5);
    end

    % Add a downward blue arrow next to the black cylinder
    arrowStart = wellCenter + [100, 0, wellHeight/2+100];  % Arrow start point above the cylinder center
    arrowEnd = arrowStart + [0, 0, -500];  % Arrow end point 500m below the start point
    quiver3(arrowStart(1), arrowStart(2), arrowStart(3), ...
            arrowEnd(1) - arrowStart(1), arrowEnd(2) - arrowStart(2), arrowEnd(3) - arrowStart(3), ...
            'b', 'LineWidth', 2, 'MaxHeadSize', 20);  % Draw a blue arrow with quiver3
    %}
    
    quiver3(100, 1950, 0, 400, 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 1.5);
    quiver3(500, 2150, 0, -400, 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 1.5);

    % Set the aspect ratio and the limits of the plot
    daspect([1 1 1]);  % Set equal data unit length on all axes
 
    % Remove axes and labels
    xlabel('');
    ylabel('');
    zlabel('');

    axis tight;  % Shrink axes limits to fit all data
    axis off;    % Remove axis lines and labels
    view(315, 30);
end
