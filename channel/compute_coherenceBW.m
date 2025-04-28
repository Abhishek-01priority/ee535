function [pdpcobw rmscobw] = compute_coherenceBW(pdp, B)

	% calculating coherence BW from pdp
	R_f = abs(fft(pdp));
	R_f = R_f / max(R_f);  % Normalize
	threshold = 0.5;
	index = find(R_f < threshold, 1, 'first');
	df = B / length(pdp);;  % frequency resolution
	pdpcobw = (index - 1) * df;
	f_axis = linspace(-B/2, B/2, length(R_f));  % Frequency axis
	figure;
	plot(f_axis/1e6, fftshift(R_f));  % in MHz
	hold on
	plot(f_axis/1e6, 0.5 * ones(size(f_axis)),'k-')
	xlabel('Frequency Offset Δf (MHz)');
	ylabel('|R_f(Δf)|');
	title(['Frequency Correlation vs Δf, B_c ≈ ' num2str(pdpcobw/1e6, '%.2f') ' MHz']);
	grid on;

	% coherence BW from RMS delay spread
	tau = (0:length(pdp)-1) * 1/B;  % Delay axis
	P = pdp / sum(pdp);       % Normalize
	tau_mean = sum(tau .* P);
	tau_rms = sqrt(sum(P .* (tau - tau_mean).^2));
	rmscobw = 1 / (5 * tau_rms);

end