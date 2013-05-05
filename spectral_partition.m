function [] = spectral_partition(A)
%SPECTRAL_PARTITION Summary of this function goes here
%   Detailed explanation goes here

L = laplacian_matrix(full(A));
[V D] = eigs(L, 2, 'sa');
disp(D(2, 2)); % second smallest eigenvector
plot(sort(V(:,2)), '.-');
[ignore p] = sort(V(:, 2));
spy(A(p, p));

end

% a1 : 0.8992
% a2 : 0.5614
% a3 : 0.9856
% a4 : 0.9116
% a5 : 0.9352

