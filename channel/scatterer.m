%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective: 	Calculate the location of scatterers such that the overall statistical properties of
%%				channel follows Nakagami + log normal distribution
%% Inputs:		TBD
%% Output:		(1)h_l: Complex channel gain from each scatterer N_s x 1 vector
%%				(2)scatterers: N_s x 2 matrix with column1->x coordinate and column2 -> y coordinate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H, h_l, scatterers, aoa, aod, t_l] = scatterer(N_s)
	%% Parameters
	d = 150;            % Tx-Rx distance in meters
	a = d/2 + 10;       % Semi-major axis
	b = 40;             % Semi-minor axis
	Nrxtx = 4; 			% number of rx = tx antennas
	lambda = 0.1; 		% wavelength
	spacing = lambda/2;		% spacing between array elements

	m = 3;            % Nakagami-m parameter
	Omega = 1;          % Average Nakagami power
	sigma_dB = 6;       % Lognormal shadowing std dev (in dB)

	% Tx and Rx positions
	Tx = [-d/2, 0];
	Rx = [d/2, 0];

	% scatter locations
	scatterers = zeros(N_s, 2);

	% angle of arrival and angle of departure for each scatterer
	aoa = zeros(N_s,1);
	aod = zeros(N_s,1);

	% rx and tx array response for a given AoA and AoD respectively
	a_r = zeros(Nrxtx, N_s);
	a_t = zeros(Nrxtx, N_s);

	% complex channel gain
	h_l = zeros(N_s, 1);

	% time delat
	t_l = zeros(N_s, 1);

	%% Main loop
    count = 0;
    while count < N_s
        % Generate scatterer within ellipse
        x = (2 * rand - 1) * a;
        y = (2 * rand - 1) * b;
        if (x^2)/(a^2) + (y^2)/(b^2) <= 1

			scatterers(count+1, :) = [x, y]; % scatterer locations

            % Nakagami fading (Rayleigh if m = 1)
            g = gamrnd(m, Omega/m);     % Nakagami power
            r_mag = sqrt(g);

            % Lognormal shadowing
            chi_dB = sigma_dB * randn;
            shadowing = 10^(chi_dB / 20);

            % random phase induced by each scatterer
            phi = 2*pi*rand;

            % Complex channel gain
            h_l(count+1) = r_mag * shadowing * exp(1j*phi);

            % angle of arrival and array response at RX
            aoa(count+1) = atan2(Rx(2) - y, Rx(1) - x); % in radians
            a_r(:, count+1) = exp(1i * 2*pi* spacing * (0:Nrxtx-1)'/lambda * sin(aoa(count+1)));

            % angle of departure and array response at TX
            aod(count+1) = atan2(y - Tx(2), x - Tx(1)); % in radians
            a_t(:, count+1) = exp(1i * 2*pi* spacing * (0:Nrxtx-1)'/lambda * sin(aod(count+1)));

            % time delay due to each scatterer
            dist = sqrt(sum((Tx - scatterers(count+1, :)).^ 2)) + ...
            	   sqrt(sum((Rx - scatterers(count+1, :)) .^ 2)); % run length calculation
            t_l(count+1) = dist / (3e8);

            % complex channel gain
            H(:,:,count+1) = h_l(count+1) * a_r(:,count+1) * a_t(:,count+1)';

            count = count + 1;
        end
    end

	% Optional: Plot locations comment/uncomment to plot this
	% figure;
	% plot(Tx(1), Tx(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
	% plot(Rx(1), Rx(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
	% scatter(scatterers(:,1), scatterers(:,2), 20, 'k', 'filled');
	% legend('Tx', 'Rx', 'Scatterers');
	% xlabel('x [m]'); ylabel('y [m]');
	% title('Scatterer locations');
	% axis equal; grid on;

	% % Plotting the AoA for each scatterer. It can be seen most of the angles are in right half 
	% % plane indicating the scatterers are in front of receiver
	% figure;
	% plot(cos(aoa), sin(aoa), 'ro', 'MarkerFaceColor', 'r')
	% hold on
	% plot(cos(aod), sin(aod), 'b^', 'MarkerFaceColor', 'b')
	% xlim([-1.2 1.2])
	% ylim([-1.2 1.2])

	% keyboard

end