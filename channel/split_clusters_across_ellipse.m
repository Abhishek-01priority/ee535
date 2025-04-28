%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% © Abhishek Manjunath 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N_cluster_per_ellipse] = split_clusters_across_ellipse(N_circles, N_clusters)

	b = 1; % decay factor
	w = zeros(1, N_circles); % weights 

	w = exp(-(1:N_circles) * b); % number of clusters should decay exponentially
	w = w / sum(w); % normalizing

	N_cluster_per_ellipse = round(w * N_clusters) + 1;

	d = N_clusters - sum(N_cluster_per_ellipse);

	[~, idx] = max(N_cluster_per_ellipse);

	N_cluster_per_ellipse(idx) += d;
end