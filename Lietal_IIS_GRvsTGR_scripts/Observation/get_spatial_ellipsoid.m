function [Lmax, Lmedian, Lmin, V, coverage] = get_spatial_ellipsoid(X, Y, Z, plt_flag)
    % Combine into matrix
    P = [X(:)'; Y(:)'; Z(:)'];
    N = size(P, 2);

    % Step 1: PCA (Eigen decomposition of covariance)
    mu = mean(P, 2);
    Sigma = cov(P');
    [V_pca, ~] = eig(Sigma);

    % Step 2: Project points onto principal axes
    P_centered = P - mu;
    P_proj = V_pca' * P_centered;

    % Step 3: Z-score filtering along each axis
    is_inlier = true(1, N);
    z_thresh = 3;  % Default Z-score threshold
    for i = 1:3
        coord = P_proj(i, :);
        z_scores = (coord - mean(coord)) / std(coord);
        is_inlier = is_inlier & (abs(z_scores) <= z_thresh);
    end

    % Step 4: Keep only inliers
    P_clean = P(:, is_inlier);
    T_clean = linspace(0, 1, numel(P_clean(2,:)));
    P_outlier = P(:, ~is_inlier);

    % Step 5: Fit ellipsoid to cleaned data
    mu_clean = mean(P_clean, 2);
    Sigma_clean = cov(P_clean');

    % Step 6: Compute Mahalanobis distances and max radius
    invSigma_clean = inv(Sigma_clean);
    d2_all = sum((P_clean - mu_clean) .* (invSigma_clean * (P_clean - mu_clean)), 1);
    d2_max = max(d2_all);  % Ensures 100% of cleaned data is enclosed

    % Step 7: Axes lengths and volume
    [Vaxes, D] = eig(Sigma_clean);
    axes_lengths = sqrt(diag(D) * d2_max);
    [Lsorted, ~] = sort(axes_lengths, 'descend');
    Lmax = Lsorted(1);
    Lmedian = Lsorted(2);
    Lmin = Lsorted(3);
    V = (4/3) * pi * Lmax * Lmedian * Lmin;

    % Step 8: Compute coverage
    coverage = sum(is_inlier) / N;

    % Step 9: Visualization
    if plt_flag == 1
        hold on; grid on;

        % Plot inliers colored by T_clean
        scatter3(P_clean(1,:), P_clean(2,:), P_clean(3,:), 30, T_clean, 'filled');
        colormap(turbo);
        caxis([0 1]);        
        
        % Plot outliers as black 'X'
        scatter3(P_outlier(1,:), P_outlier(2,:), P_outlier(3,:), 10, 'k', 'x', 'LineWidth', 1.5);
        

        % Draw ellipsoid in gray
        drawEllipsoid(mu_clean, Sigma_clean, d2_max, [0.5 0.5 0.5]);

        axis equal;
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        view(3); camlight; lighting gouraud;
        set(gca,'fontsize',16);
        
    end
end

function drawEllipsoid(mu, Sigma, scale, color)
    [V, D] = eig(Sigma);
    [x, y, z] = ellipsoid(0, 0, 0, 1, 1, 1, 40);
    XYZ = [x(:) y(:) z(:)]';
    A = V * sqrt(D * scale);
    ellip = A * XYZ + mu;
    x = reshape(ellip(1,:), size(x));
    y = reshape(ellip(2,:), size(y));
    z = reshape(ellip(3,:), size(z));
    surf(x, y, z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', color);
end

