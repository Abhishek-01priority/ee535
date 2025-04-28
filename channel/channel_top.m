%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective: 	Top level code for running channel simulation. This is not a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

geoparams = struct();
geoparams.c = 50;        % distance of TX and RX from origin
geoparams.R_circle = [75 150 200]; % radius of circle of coverage

B = 20e6; % channel bandwidth
Ts = 1/B; % sampling interval

nfft = 1024;
N_clusters = 80;
N_s = 80; %25;          % Scatterers per realization
N_realizations = 1000;

N_taps = floor( sqrt( geoparams.c^2 + geoparams.R_circle(end)^2 ) * 2 / 3e8 * B ) + 10; % 9 extra taps for safety

h_l = aoa = aod = zeros(N_s, N_realizations);
Hf = zeros(4,4,nfft,N_realizations);
cir = zeros(4,4, N_taps, N_realizations); % CIR be slighly larger than max observed delay and for every realization
envelopes = zeros(N_realizations,1);

%% Main loop
tic()
for r = 1:N_realizations

	% caculate channel params
    % [H, h_l, scatterers, aoa, aod, t_l] = scatterer(N_s);
    % [H, h_l, scatterers, aoa, aod, t_l] = concentric_ellipse(N_s);
    [H, h_l(:, r), scatterers, aoa(:, r), aod(:, r), t_l] = circular_clusters(N_clusters, N_s / N_clusters, geoparams);

    % calculate channel impulse response by aligning complex channel gains with respective time delays
    % if idx > N_taps, error will be hit
    idx = floor(t_l / Ts) + 1;
    for i = 1 : N_s
    	cir(:, :, idx(i), r) = cir(:, :, idx(i), r) + H(:,:,i);
    end
    envelopes(r) = abs(sum( squeeze( cir(1, 1, :, r) ) )); % fetch first tx to first rx CIR for a given channel realization

    % Calculating the frequency response
    % Hf(:,:,:,r) = fft(cir(:,:,:,r), nfft, 3);

    % Calculating the frequency response
    t1 = permute(cir(:,:,:,r), [3 1 2]);
    t2 = reshape(t1, N_taps, []);
    t3 = fft(t2, nfft, 1);
    t4 = reshape(t3, nfft, 4, 4);
    t5 = permute(t4, [2 3 1]);
    Hf(:, :, :, r) = t5;

    pause(0.01);
    progress_bar(r, N_realizations);

    % keyboard

end
fprintf("\n");
toc()

% %% Plot histogram
% figure;
% [counts, edges] = hist(envelopes, 50);      % 50 bins
% bin_width = edges(2) - edges(1);
% pdf_vals = counts / (sum(counts) * bin_width);  % Normalize to PDF

% plot(edges, pdf_vals, 'k', 'LineWidth', 1.5, 'DisplayName', 'hist of channel');
% hold on;

% % Fit Nakagami distribution
% pd_naka = fitdist(envelopes, 'nakagami');
% x = linspace(min(envelopes), max(envelopes), 500);
% y_naka = pdf(pd_naka, x);
% plot(x, y_naka, 'r', 'LineWidth', 1.8, 'DisplayName', 'Fitted Nakagami');

% % Fit lognormal distribution
% pd_logn = fitdist(envelopes, 'lognormal');
% y_logn = pdf(pd_logn, x);
% plot(x, y_logn, 'g-', 'LineWidth', 1.8, 'DisplayName', 'Fitted Lognormal');

% xlabel('|h|'); ylabel('PDF');
% title('Nakagami-Lognormal Composite Fading (Envelope)');
% legend('Location','northeast');
% grid on;

% figure;
% plot(10*log10(envelopes))
% xlabel('Channel realizations')
% ylabel('channel gain in dBm')