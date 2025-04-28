function [R_rx, R_tx] = compute_correlation_matrices(H_tap)
  [Nr, Nt, N_taps, N_real] = size(H_tap);

  % Reshape H_tap to [Nr x Nt x (N_taps * N_real)]
  H_all = reshape(H_tap, Nr, Nt, []);

  % Initialize
  R_rx = zeros(Nr, Nr);
  R_tx = zeros(Nt, Nt);

  % Now vectorized summation
  for idx = 1:size(H_all, 3)
    Hn = H_all(:,:,idx);
    R_rx += Hn * Hn';
    R_tx += Hn' * Hn;
  end

  % Normalize
  R_rx /= size(H_all,3);  % Same as (N_taps * N_real)
  R_tx /= size(H_all,3);

  % Plot
  figure;
  subplot(1,2,1); imagesc(abs(R_rx)); colorbar;
  title('Receive Correlation Matrix');
  subplot(1,2,2); imagesc(abs(R_tx)); colorbar;
  title('Transmit Correlation Matrix');
end
