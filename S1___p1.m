% Load the necessary data for camera intrinsic matrix and point correspondences
load('data/intrinsics.mat'); % Contains the camera intrinsic matrix K
load('data/some_corresp.mat'); % Contains matched points between two images

% Define image paths for the stereo images
img1Path = 'data/im1.png'; % Replace with actual path
img2Path = 'data/im2.png'; % Replace with actual path

% Load the stereo images
img1 = imread(img1Path);
img2 = imread(img2Path);

% Set M as the max dimension from both images, which is used for normalizing points
M = max([size(img1, 1), size(img1, 2), size(img2, 1), size(img2, 2)]);

% Call the eight_point function with loaded correspondences to compute the Fundamental matrix
F = eight_point(pts1, pts2, M);


% Randomly select a few points from pts1
num_points_to_visualize = 8;
random_indices = randperm(size(pts1, 1), num_points_to_visualize);
selected_pts1 = pts1(random_indices, :);

% Compute corresponding epipolar lines in the second image
epiLines = epipolarLine(F', selected_pts1);

% Compute the intersection points of the epipolar lines with the image borders
points = lineToBorderPoints(epiLines, size(img2));


% Plot the epipolar lines on the second image
figure;
imshow(img2);
title('Epipolar Lines in Second Image');
hold on;

% Plot the epipolar lines for selected points
for i = 1:size(points, 1)
    plot(points(i, [1,3])', points(i, [2,4])', 'g'); % Plotting line between border points
end

% Scatter plot of the points in the first image corresponding to the epipolar lines
scatter(selected_pts1(:,1), selected_pts1(:,2), 'ro');
hold off;



% Function to normalize points in preparation for the Eight-Point Algorithm
function [normalized_points, T] = normalize_points(points, M)
    % Convert points to double precision for arithmetic operations
    points = double(points);
    
    % Compute the centroid of the points
    centroid = mean(points, 1);
    
    % Translate points so the centroid is at the origin
    translated_points = points - centroid;
    
    % Scale points so that the average distance from the origin is sqrt(2)
    scaling_factor = sqrt(2) / mean(sqrt(sum(translated_points.^2, 2)));
    
    % Create transformation matrix to normalize points
    T = [scaling_factor, 0, -scaling_factor * centroid(1);
         0, scaling_factor, -scaling_factor * centroid(2);
         0, 0, 1];
    
    % Apply transformation to points
    normalized_points = (T * [points, ones(size(points, 1), 1)]')';
    normalized_points = normalized_points(:, 1:2);
end

% Function to construct matrix A used in the Eight-Point Algorithm
function A = construct_matrix_A(pts1, pts2)
    % Number of points correspondences
    N = size(pts1, 1);
    
    % Construct matrix A as defined in the Eight-Point Algorithm
    A = [pts2(:, 1) .* pts1(:, 1), pts2(:, 1) .* pts1(:, 2), pts2(:, 1), ...
         pts2(:, 2) .* pts1(:, 1), pts2(:, 2) .* pts1(:, 2), pts2(:, 2), ...
         pts1(:, 1), pts1(:, 2), ones(N, 1)];
end

% Main function to implement the Eight-Point Algorithm
function F = eight_point(pts1, pts2, M)
    % Normalize the points
    [pts1_norm, T1] = normalize_points(pts1, M);
    [pts2_norm, T2] = normalize_points(pts2, M);
    
    % Construct the matrix A using normalized points
    A = construct_matrix_A(pts1_norm, pts2_norm);

    % Compute the fundamental matrix F using SVD on matrix A
    [~, ~, V] = svd(A);
    F = reshape(V(:, end), 3, 3)';
    
    % Enforce the rank-2 constraint by setting the smallest singular value to 0
    [Uf, Sf, Vf] = svd(F);
    Sf(3, 3) = 0;
    F = Uf * Sf * Vf';
    
    % Denormalize the fundamental matrix using the transformation matrices
    F = T2' * F * T1;
end

