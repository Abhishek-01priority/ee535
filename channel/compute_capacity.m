function capacity_curve = compute_capacity(H_tap, SNR_dB)
  [Nr, Nt, N_taps, N_real] = size(H_tap);
  SNR_lin = 10.^(SNR_dB/10);

  % Flatten taps and realizations together
  H_all = reshape(H_tap, Nr, Nt, []);  % [Nr x Nt x (N_taps*N_real)]
  N_matrices = size(H_all, 3);

  capacity_curve = zeros(length(SNR_dB), 1);

  for k = 1:length(SNR_dB)
    snr_factor = SNR_lin(k)/Nt;
    cap_sum = 0;

    for idx = 1:N_matrices
      Hn = H_all(:,:,idx);
      cap_sum += real(log2(det(eye(Nr) + snr_factor * (Hn * Hn'))));
    end

    capacity_curve(k) = cap_sum / N_matrices;
  end

  % Plot
  figure;
  plot(SNR_dB, capacity_curve, '-o');
  xlabel('SNR (dB)'); ylabel('Average Capacity (bps/Hz)');
  title('MIMO Capacity vs SNR');
  grid on;
end
