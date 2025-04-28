%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% © Abhishek Manjunath 2025
%% inputs:  aoa, aod, alpha -> N_s x N_realizations
%% outputs: sigma_aoa, sigma_aod -> 1 x N_realizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma_aoa, sigma_aod] = compute_angular_spread(aoa, aod, alpha)
  
  weights = abs(alpha).^2;
  aoa_mean = sum(weights .* aoa) / sum(weights, 1);
  aod_mean = sum(weights .* aod) / sum(weights, 1);

  sigma_aoa = sqrt(sum(weights .* (aoa - aoa_mean).^2, 1) / sum(weights, 1));
  sigma_aod = sqrt(sum(weights .* (aod - aod_mean).^2, 1) / sum(weights, 1));

end