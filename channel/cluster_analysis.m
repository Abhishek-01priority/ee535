function cluster_analysis(phi_r, phi_t, alpha)
  % Inputs:
  % phi_r, phi_t : AoA/AoD (L x 1 vectors), radians
  % alpha        : complex path gains (L x 1 vector)
  L = length(alpha);

  % Convert angles to degrees for intuitive interpretation
  X = [rad2deg(phi_t), rad2deg(phi_r)];  % [AoD, AoA]

  % Choose number of clusters (e.g., 3–5 initially, adjust as needed)
  K = 4;

  % K-means clustering
  [idx, centers] = kmeans(X, K, 'Replicates', 10);

  % Compute total power per cluster
  cluster_power = zeros(1,K);
  for k = 1:K
    cluster_power(k) = sum(abs(alpha(idx == k)).^2);
  end

  % Sort clusters by power
  [cluster_power, sort_idx] = sort(cluster_power, 'descend');
  centers = centers(sort_idx,:);

  % Display cluster analysis
  fprintf('Cluster Analysis:\n');
  for k = 1:K
    fprintf('Cluster %d: AoD = %.1f°, AoA = %.1f°, Power fraction = %.2f%%\n', ...
        k, centers(k,1), centers(k,2), 100*cluster_power(k)/sum(cluster_power));
  end

  % Visualization
  figure;
  scatter(X(:,1), X(:,2), 50, idx, 'filled');
  hold on;
  scatter(centers(:,1), centers(:,2), 100, 'kx', 'linewidth', 2);
  xlabel('AoD (degrees)'); ylabel('AoA (degrees)');
  title('AoA/AoD Cluster Analysis');
  grid on;
  colorbar;

end