function [aoa_spread, aod_spread] = aps(phi_r, phi_t, alpha)

	% Parameters
	N_bins = 36;  % 10° angular resolution
	edges = linspace(-pi, pi, N_bins+1);
	bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
	alpha_power = abs(alpha).^2;

	% Initialize heatmap matrix
	aps_map = zeros(N_bins, N_bins);  % rows: AoA bins, cols: AoD bins

	for l = 1:length(alpha)
	  i = find(phi_r(l) >= edges(1:end-1) & phi_r(l) < edges(2:end));
	  j = find(phi_t(l) >= edges(1:end-1) & phi_t(l) < edges(2:end));
	  if ~isempty(i) && ~isempty(j)
	    aps_map(i, j) += alpha_power(l);
	  end
	end

	% Normalize
	aps_map /= sum(aps_map(:));

	% % Plot heatmap
	% figure;
	% imagesc(rad2deg(bin_centers), rad2deg(bin_centers), aps_map);
	% xlabel('AoD (degrees)');
	% ylabel('AoA (degrees)');
	% title('Joint Angular Power Spectrum (Heatmap)');
	% colorbar;
	% axis xy;

	% Array Parameters
	Nr = 4; Nt = 4; lambda_c = 1; d = lambda_c/2;

	% RMS angular spread function
	function spread = rms_ang_spread(phi, power)
	  mean_angle = angle(sum(power .* exp(1j*phi)));
	  angular_deviation = angle(exp(1j*(phi - mean_angle)));
	  spread = sqrt(sum(power .* angular_deviation.^2) / sum(power));
	  spread = rad2deg(spread);
	end

	% Compute angular spreads
	aoa_spread = rms_ang_spread(phi_r, alpha_power);
	aod_spread = rms_ang_spread(phi_t, alpha_power);

	% Array response vector function
	function a = array_resp(N, angle, d, lambda_c)
	  n = (0:N-1)';
	  a = exp(1j*2*pi*d*n*sin(angle)/lambda_c) / sqrt(N);
	end

	% Initialize spatial correlation matrices
	R_rx = zeros(Nr,Nr);
	R_tx = zeros(Nt,Nt);

	for l = 1:length(alpha)
	  ar = array_resp(Nr, phi_r(l), d, lambda_c);
	  at = array_resp(Nt, phi_t(l), d, lambda_c);
	  R_rx += alpha_power(l) * (ar * ar');
	  R_tx += alpha_power(l) * (at * at');
	end

	% % Normalize correlation matrices
	% R_rx /= trace(R_rx);
	% R_tx /= trace(R_tx);

	% % Plot spatial correlation matrices
	% figure;
	% subplot(1,2,1);
	% imagesc(abs(R_rx)); colorbar;
	% title(sprintf('Rx Spatial Correlation\nAoA spread=%.2f°', aoa_spread));
	% xlabel('Rx Antenna Index'); ylabel('Rx Antenna Index');

	% subplot(1,2,2);
	% imagesc(abs(R_tx)); colorbar;
	% title(sprintf('Tx Spatial Correlation\nAoD spread=%.2f°', aod_spread));
	% xlabel('Tx Antenna Index'); ylabel('Tx Antenna Index');


end