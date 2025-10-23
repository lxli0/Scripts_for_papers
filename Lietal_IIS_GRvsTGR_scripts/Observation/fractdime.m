function df = fractdime(x, y, z, plt_flag)
    n = length(x);
    
    % Compute upper triangular pairwise distances efficiently
    R = pdist([x(:), y(:), z(:)]);  % 3D distances

    % Logarithmic binning using linspace for radius
    hk = 1000;
    rmin = min(R);
    rmax = max(R);
    rc_edges = linspace(rmin, rmax, hk);
    rc = (rc_edges(1:end-1) + rc_edges(2:end)) / 2;

    % Histogram of distances
    counts = histcounts(R, rc_edges);

    % Normalized correlation sum
    norm_factor = n * (n - 1) / 2;
    c = 2 * cumsum(counts) / norm_factor;

    % Log-log plot for correlation dimension
    lgrc = log10(rc);
    lgc = log10(c);

    % Fit in appropriate scaling region
    mask = c >= 5*min(c) & c <= 0.2;
    p = polyfit(lgrc(mask), lgc(mask), 1);
    df = p(1);

    % Optional plotting
    if plt_flag == 1
        scatter(lgrc, lgc, 30, 1/255 * [37 122 182], 'filled'); hold on;
        plot(lgrc(mask), polyval(p, lgrc(mask)), 'color', 1/255 * [252 132 13], 'linewidth', 2.5);
        xlabel('$\log r$', 'Interpreter', 'latex');
        ylabel('$\log C$', 'Interpreter', 'latex');
        box on; grid on; grid minor;
        set(gca, 'fontsize', 16);
        textStr = sprintf('$d_f$: %.2f', df);
        xLimits = xlim;
        yLimits = ylim;
        text(xLimits(1) + 0.05 * diff(xLimits), ...
             yLimits(2) - 0.05 * diff(yLimits), ...
             textStr, 'Interpreter', 'latex', 'FontSize', 16);
    end
end
