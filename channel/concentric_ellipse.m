function [H, alpha, scatterers, aoa, aod, tau] = concentric_ellipse(N_s)

	% Parameters
	c = 150; tx = [-c, 0]; rx = [c, 0];
	a_list = [180, 200, 220];
	area_list = pi * a_list .* sqrt(a_list.^2 - c^2);
	total_area = sum(area_list);
	N_per_ellipse = round(N_s * (area_list / total_area));
	% N_per_ellipse = N_s/length(a_list); % same number of scatterers per ellipse
	lambda = 1; d = lambda/2;
	Nt = 4; Nr = 4;
	c_light = 3e8;
	m_nak = 1.5; mu_db = 0; sigma_db = 4;

	% Initialization
	L = length(a_list) * N_s;
	scatterers = zeros(L, 2);
	aod = zeros(L, 1); aoa = zeros(L, 1);
	tau = zeros(L, 1); alpha = zeros(L, 1);
	At = zeros(Nt, L); Ar = zeros(Nr, L);

	idx = 1;
	for i = 1:length(a_list)
	  a = a_list(i);
	  b = sqrt(a^2 - c^2);

	  for k = 1:N_per_ellipse(i)
	    while true
	      r = sqrt(rand()); theta = 2*pi*rand();
	      x = r*a*cos(theta); y = r*b*sin(theta);
	      if (x^2)/(a^2) + (y^2)/(b^2) <= 1
	        break;
	      end
	    end

	    scatterers(idx,:) = [x, y];

	    % Angles
	    vec_tx = [x, y] - tx;
	    vec_rx = rx - [x, y];
	    aod(idx) = atan2(vec_tx(2), vec_tx(1));
	    aoa(idx) = atan2(vec_rx(2), vec_rx(1));

	    % Delay
	    d_tx = norm([x, y] - tx);
	    d_rx = norm([x, y] - rx);
	    tau(idx) = (d_tx + d_rx) / c_light;

	    % Gain
	    nak = sqrt(gamrnd(m_nak, 1/m_nak));
	    logn = 10^(normrnd(mu_db, sigma_db)/20);
	    alpha(idx) = nak * logn * exp(1j * 2*pi*rand());

	    % Array response
	    n_tx = (0:Nt-1)'; n_rx = (0:Nr-1)';
	    At(:,idx) = exp(1j * 2*pi*d * n_tx * sin(aod(idx))) / sqrt(Nt);
	    Ar(:,idx) = exp(1j * 2*pi*d * n_rx * sin(aoa(idx))) / sqrt(Nr);

	    H(:, :, idx) = alpha(idx) * At(:, idx) * Ar(:, idx)';

	    idx += 1;
	  end
	end

	% % Plot: Geometry
	% theta = linspace(0, 2*pi, 300);
	% figure; hold on; axis equal; grid on;
	% colors = {'b', 'g', 'r'};
	% for i = 1:length(a_list)
	%   a = a_list(i); b = sqrt(a^2 - c^2);
	%   x = a * cos(theta); y = b * sin(theta);
	%   plot(x, y, colors{i});
	% end
	% scatter(scatterers(:,1), scatterers(:,2), 20, 'k', 'filled');
	% plot(tx(1), tx(2), 'bs', 'MarkerFaceColor', 'b');
	% plot(rx(1), rx(2), 'r^', 'MarkerFaceColor', 'r');
	% text(tx(1)-0.5, tx(2)-0.5, 'Tx');
	% text(rx(1)+0.5, rx(2)-0.5, 'Rx');
	% title('Scatterer Placement Inside Multiple Ellipses');
	% xlabel('x-axis'); ylabel('y-axis');

	% keyboard

end