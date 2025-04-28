%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% © Abhishek Manjunath 2025
%% Objective:	a) Power Delay Profile
%% 				b) Coherence Bandwidth
%% 				c) Rank
%%				d) Eigen Spectrum
%% 				e) angular spread
%% 				f) correlation matrix for rx and tx antenna
%% 				g) Capacity
%% Inputs:		a) Channel impulse response N x N x N_Taps x N_realizations
%% 				b) Frequeny response of channel H(f)
%% 				c) Complex channel gain
%% 				d) Sampling time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analysis_of_channel(cir, H_f, aoa, aod, h_l, Ts)

	%% a) Power Delay profile
	[pdp, delay_us] = compute_pdp(cir, Ts);
	fprintf("Maximum delay = %0.2f us\n", max(delay_us(find(delay_us ~= 0))));
	% ensuring PDP stays a column vector
	if (size(pdp, 2) == 1)
		pdp = reshape(pdp, 1 ,[]);
	end

  	%% b) Coherence Bandwidth
  	[cobw, rms_cobw] = compute_coherenceBW(pdp, 1/Ts);
  	fprintf("Coherence BW from pdp = %0.2f MHz, from RMS delay spread = %0.2f MHz\n", cobw/1e6, rms_cobw/1e6);

  	%% c) Rank d) Eigen values
  	[avg_rank_per_subcarrier, condition_no, eigenvals] = compute_rank(H_f);

  	%% e) angular spread
  	[sigma_aoa, sigma_aod] = compute_angular_spread(aoa, aod, h_l);
	fprintf('RMS AoA Spread = %.2f degrees\n', rad2deg(sigma_aoa));
	fprintf('RMS AoD Spread = %.2f degrees\n', rad2deg(sigma_aod));

	%% f) Correlation matrices
	[R_rx, R_tx] = compute_correlation_matrices();

	%% g) capacity

end