function [r_avg, k_number, eigenvals] = compute_rank(H_f)
	
	step_size = 4; % added this because SVD is inherently slow for 4x4 and it was taking too long
	k_number = zeros(size(H_f, 3), size(H_f, 4));
	r_avg = zeros(1, size(H_f, 3));
	eigenvals = zeros(size(H_f, 1), size(H_f, 3));

	% iterate through different channel realizations
	for r = 1 : step_size : size(H_f, 4)

		% go over all the subcarrier to compute rank for each subcarrier
		for n = 1 : size(H_f, 3)

			H_per_subcarrier = H_f(:,:,n,r);

			% rank estimation
			s = svd(H_per_subcarrier);
			r_avg(n) += length(find(s > 0.01 * max(s))); % taking 1% as threshold

			% condition number per subcarrier
			k_number(n,r) = max(s) / min(s);

			% eigenvalues
			eigenvals(:, n) += s.^2;
		end
	end

	r_avg = step_size/size(H_f, 4) * r_avg;
	k_number(k_number == 0) = [];
	k_number = log10(k_number);
	eigenvals = eigenvals / size(H_f, 4);

	% figure;
	% stem(1:size(H_f,3), r_avg, 'filled')
	% xlabel('Subcarrier Index');
  	% ylabel('Average Rank');
  	% title('Average Rank vs Subcarrier');
  	% grid on;

end