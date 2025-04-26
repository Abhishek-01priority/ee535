%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% © Abhishek Manjunath 2025
%% Learning:  1) Sometimes the cluster centers overlap but it should be completely fine because 
%%               they will have same angular spread and phyiscally can be interpreted as same cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H, alpha, scatterers, aoa, aod, tau] = circular_clusters(N_clusters, N_scatters_per_cluster)

  % Geometry Setup
  c = 150;                   % Half-distance between Tx and Rx
  tx = [-c, 0];
  rx = [c, 0];
  R = 200;                   % Radius of overall circle
  cluster_spread = 10;        % Std-dev spread within cluster (meters)

  lambda = 1;
  d = lambda/2;
  Nt = 4; Nr = 4;
  c_light = 3e8;

  m_nak = 1.5; mu_db = 0; sigma_db = 4;

  total_L = N_clusters * N_scatters_per_cluster;
  scatterers = zeros(total_L, 2);
  cluster_centers = zeros(N_clusters, 2);
  aoa = zeros(total_L, 1);
  aod = zeros(total_L, 1);
  tau = zeros(total_L, 1);
  alpha = zeros(total_L, 1);
  At = zeros(Nt, total_L);
  Ar = zeros(Nr, total_L);

  % Step 1: Generate cluster centers randomly inside a big circle
  for c_idx = 1:N_clusters
    while true
      x_center = (2*rand()-1)*R;
      y_center = (2*rand()-1)*R;
      if x_center^2 + y_center^2 <= R^2
        cluster_centers(c_idx,:) = [x_center, y_center];
        break;
      end
    end
  end

  % Step 2: Generate scatterers around each cluster center
  idx = 1;
  for c_idx = 1:N_clusters
    center = cluster_centers(c_idx,:);
    for k = 1:N_scatters_per_cluster
      offset = cluster_spread * randn(1,2);
      pos = center + offset;
      scatterers(idx,:) = pos;

      % Angles
      vec_tx = pos - tx;
      vec_rx = rx - pos;
      aod(idx) = atan2(vec_tx(2), vec_tx(1));
      aoa(idx) = atan2(vec_rx(2), vec_rx(1));

      % Delay
      d_tx = norm(pos - tx);
      d_rx = norm(pos - rx);
      tau(idx) = (d_tx + d_rx) / c_light;

      % Gain
      nak = sqrt(gamrnd(m_nak, 1/m_nak));
      logn = 10^(normrnd(mu_db, sigma_db)/20);
      alpha(idx) = nak * logn * exp(1j * 2*pi*rand());

      % Array responses
      n_tx = (0:Nt-1)';
      n_rx = (0:Nr-1)';
      At(:,idx) = exp(1j * 2*pi*d * n_tx * sin(aod(idx))) / sqrt(Nt);
      Ar(:,idx) = exp(1j * 2*pi*d * n_rx * sin(aoa(idx))) / sqrt(Nr);

      % channel matrix
      H(:, :, idx) = alpha(idx) * (Ar(:,idx) * At(:,idx)');

      idx += 1;
    end
  end


  % % Plot Geometry
  % theta = linspace(0, 2*pi, 300);
  % figure; hold on; axis equal; grid on;
  % plot(R*cos(theta), R*sin(theta), 'r--', 'DisplayName', 'Boundary Circle');

  % scatter(cluster_centers(:,1), cluster_centers(:,2), 100, 'g', 'filled', 'DisplayName', 'Cluster Centers');
  % scatter(scatterers(:,1), scatterers(:,2), 15, 'k', 'filled', 'DisplayName', 'Scatterers');
  % plot(tx(1), tx(2), 'bs', 'MarkerFaceColor', 'b', 'DisplayName', 'Tx');
  % plot(rx(1), rx(2), 'r^', 'MarkerFaceColor', 'r', 'DisplayName', 'Rx');
  % text(tx(1)-10, tx(2)-10, 'Tx');
  % text(rx(1)+10, rx(2)-10, 'Rx');
  % title('Scatterers Grouped in Clusters Inside Circle');
  % xlabel('x-axis (m)'); ylabel('y-axis (m)');
  % legend('show');
end
