%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% © Abhishek Manjunath 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pdp, delay_us] = compute_pdp(cir, Ts)
	% for each tap and realization, calculate the total power of N x N matrix
	[~, ~, N_taps, ~] = size(cir);
	cirsum = sum( sum( abs(cir).^2 , 1), 2); % first sum across N_r and then across N_t, result will
											 % be 1 x 1 x N_taps X N_realizations
	pdp = mean( squeeze(cirsum), 2); % take mean across the channel realizations. squeeze() will
									 % remove the first and second dimension and result will be 
									 % N_taps x N_realizations
	delay_us = [0 : N_taps-1] * Ts * 1e6; % finding delays in nanoseconds. will be use full in calculating
										  % symbol time
	figure;
	stem(delay_us, pdp, 'filled')
	xlabel('Delay (us)'); ylabel('Power');
  	title('Average Power Delay Profile');
  	grid on;

end