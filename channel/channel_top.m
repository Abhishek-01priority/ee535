%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective: 	Top level code for running channel simulation. This is not a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

B = 20e6; % channel bandwidth
Ts = 1/B; % sampling interval

N_s = 60; %25;          % Scatterers per realization
N_realizations = 1000;

aoa_spread = aod_spread = zeros(1, N_realizations); % rms spread
R_H = zeros(4,4); % correlation matrix
ns1 = zeros(1, N_realizations);
cir = zeros(4,4, N_s, N_realizations); % CIR be slighly larger than max observed delay and for every realization
envelopes = zeros(N_realizations,1);

%% Main loop
for r = 1:N_realizations

	% caculate channel params
    % [H, h_l, scatterers, aoa, aod, t_l] = scatterer(N_s);
    [H, h_l, scatterers, aoa, aod, t_l] = concentric_ellipse(N_s);

    % calculate channel impulse response by aligning complex channel gains with respective time delays
    idx = floor(t_l / Ts) + 1;
    for i = 1 : N_s
    	cir(:, :, idx(i), r) = cir(:, :, idx(i), r) + H(:,:,i);
    end
    envelopes(r) = abs(sum( squeeze( cir(1, 1, :, r) ) )); % fetch first tx to first rx CIR for a given channel realization

    %% Instantaneous Rank
    % find the sum of cir across all taps for each realization.
    % taps are in 3rd dimension and after sum it becomes 4x4x1x1000 dim.
    Hsum = sum(cir(:,:,:,r), 3);
    s = svd(Hsum);
    ns1(r) = sum(s > 0.1 * max(s)); % number of streams

    %% Effective Rank E[HH']
    R_H = R_H + Hsum * Hsum';

    %% Instantaneous Eigen spread
    % to identify which spatial streams
    es1(:, r) = sort(abs(eig(Hsum * Hsum')), 'descend');

    %% Angular Power spectrum
    % [aoa_spread(r), aod_spread(r)] = aps(aoa, aod, h_l);

    %% cluster analysis
    % cluster_analysis(aoa, aod, h_l);

    pause(0.01);
    progress_bar(r, N_realizations);

    % keyboard

end
fprintf("\n");

%% Effective Rank Contd.
eigenvalues = eig(R_H / N_realizations);
ns2 = sum(eigenvalues > 0.1 * max(eigenvalues)); % number of streams effectively
figure;
hist(ns1);
xlabel("Rank")
ylabel("Frequency")
title(["Instantaneous Rank over channel use. Rank over time E[HH'] = " num2str(ns2, '%d')])

%% Average Eigen spread
evag = 1/1000 * sum(es1, 2);
stem(evag)
ylabel("Eigen values")
xlabel("eigen indices")
title("Average eigen values across all channel realizations")

%% Power delay profile calculation. sum up the power in individual taps across different channel use
cirsum = sum(abs(cir).^2, 4);
data = squeeze(cirsum(1, 1, :))';
pdp = 1 / N_realizations * data;
figure;
stem(pdp)
xlabel("Channel taps")
ylabel("Magnitude (linear)")
title("Power Delay Profile")

%% coherence bandwidth
R_f = abs(fft(pdp));
R_f = R_f / max(R_f);  % Normalize
threshold = 0.5;
index = find(R_f < threshold, 1, 'first');
K = length(pdp);
df = B / K;  % frequency resolution
B_c = (index - 1) * df;
f_axis = linspace(-B/2, B/2, length(R_f));  % Frequency axis
figure;
plot(f_axis/1e6, fftshift(R_f));  % in MHz
hold on
plot(f_axis/1e6, 0.5 * ones(size(f_axis)),'k-')
xlabel('Frequency Offset Δf (MHz)');
ylabel('|R_f(Δf)|');
title(['Frequency Correlation vs Δf, B_c ≈ ' num2str(B_c/1e6, '%.2f') ' MHz']);
grid on;

%% Coherence bandwidth from RMS delay spread
tau = (0:length(pdp)-1) * Ts;  % Delay axis
P = pdp / sum(pdp);       % Normalize
tau_mean = sum(tau .* P);
tau_rms = sqrt(sum(P .* (tau - tau_mean).^2));
Bc_rms = 1 / (5 * tau_rms);


%% Plot histogram
figure;
[counts, edges] = hist(envelopes, 50);      % 50 bins
bin_width = edges(2) - edges(1);
pdf_vals = counts / (sum(counts) * bin_width);  % Normalize to PDF

plot(edges, pdf_vals, 'k', 'LineWidth', 1.5, 'DisplayName', 'hist of channel');
hold on;

% Fit Nakagami distribution
pd_naka = fitdist(envelopes, 'nakagami');
x = linspace(min(envelopes), max(envelopes), 500);
y_naka = pdf(pd_naka, x);
plot(x, y_naka, 'r', 'LineWidth', 1.8, 'DisplayName', 'Fitted Nakagami');

% Fit lognormal distribution
pd_logn = fitdist(envelopes, 'lognormal');
y_logn = pdf(pd_logn, x);
plot(x, y_logn, 'g-', 'LineWidth', 1.8, 'DisplayName', 'Fitted Lognormal');

xlabel('|h|'); ylabel('PDF');
title('Nakagami-Lognormal Composite Fading (Envelope)');
legend('Location','northeast');
grid on;

figure;
plot(10*log10(envelopes))
xlabel('Channel realizations')
ylabel('channel gain in dBm')