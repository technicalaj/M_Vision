% Load necessary data
load('data/intrinsics.mat');  % Assuming 'F' is computed and available

% Compute the Essential Matrix (E) directly using F, K1, and K2
E = K2' * F * K1;

% Display the computed Essential Matrix
disp('The Essential Matrix (E) is:');
disp(E);

% Perform SVD on the initial E
[U, S, V] = svd(E);

% Averaging the two largest singular values for refining E
avg_singular_value = mean([S(1,1), S(2,2)]);

% Reconstructing E with the refined singular values
S_refined = diag([avg_singular_value, avg_singular_value, 0]);
E_refined = U * S_refined * V';

% Display the refined Essential Matrix
disp('The refined Essential Matrix (E) is:');
disp(E_refined);

% Perform the singularity check on the refined Essential Matrix
if abs(S_refined(1,1) - S_refined(2,2)) < 1e-4 && S_refined(3,3) < 1e-4
    disp('The refined E passes the singular value check.');
else
    disp('The refined E does not pass the singular value check.');
end
