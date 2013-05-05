%% Import the data.
M = importdata('nodes.txt', '	', 0);
save 'facebook.mat' M;

%% Load data and create the adjancency matrix.
load 'facebook.mat';

I = M(:, 1);
J = M(:, 2);
N = max(max(I), max(J));
A = sparse(N, N);
for i=1:length(I)
    A(I(i), J(i)) = 1;
    A(J(i), I(i)) = 1;
    disp(i/length(I));
end

spy(A);

save 'facebook_adj.mat' A;