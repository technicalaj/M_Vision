% --- Data Loading Section ---
load('data/intrinsics.mat'); % Load camera intrinsic matrix K
load('data/some_corresp.mat'); % Load matched points between two images

img1Path = 'data/im1.png'; % Define image path for the first image
img2Path = 'data/im2.png'; % Define image path for the second image

img1 = imread(img1Path); % Load the first image
img2 = imread(img2Path); % Load the second image

% Compute the Fundamental Matrix (Assuming you have the eight_point function)
M = max([size(img1, 1), size(img1, 2), size(img2, 1), size(img2, 2)]);
F = eight_point(pts1, pts2, M);


% --- Find Correspondences ---
pts2 = epipolar_correspondences(img1, img2, F, pts1);

% --- Visualization Section ---
% Plotting the points in the first image
subplot(1, 2, 1);
imshow(img1); hold on;
plot(pts1(:, 1), pts1(:, 2), 'go', 'MarkerSize', 10, 'LineWidth', 1.5);
title('Points in First Image');
hold off;

% Plotting the corresponding points in the second image
subplot(1, 2, 2);
imshow(img2); hold on;
plot(pts2(:, 1), pts2(:, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
title('Corresponding Points in Second Image');
hold off;

% --- Epipolar Correspondences Function Definition ---
function pts2 = epipolar_correspondences(img1, img2, F, pts1)
    % Preprocess images
    img1 = im2double(img1);
    img2 = im2double(img2);

    % Initialization
    pts2 = zeros(size(pts1));
    windowSize = 5; % Example window size for SSD comparison
    halfWindowSize = floor(windowSize / 2);
    [rows, cols] = size(img2(:,:,1));
    
    for i = 1:size(pts1,1)

        
        % Convert pt1 to double to ensure compatibility for matrix multiplication
        pt1 = double([pts1(i, :) 1]');

        % Ensure F is also of type double; this is usually the case, but if F is loaded from a file or computed in a way that might lead to it being an integer, explicitly convert it as well
        F = double(F);

        % Compute the corresponding epipolar line in img2
        epLine = F * pt1;

        
        % Create search lines along the epipolar line in img2
        X2 = 1:cols;
        Y2 = round((-epLine(3) - epLine(1)*X2) / epLine(2));
        
        % Initialize minimum SSD and its index
        minSSD = inf;
        bestIndex = 1;
        
        for j = 1:length(X2)
            x = X2(j);
            y = Y2(j);
            
            % Ensure the window is fully within the image bounds
            if x <= halfWindowSize || x >= cols-halfWindowSize || ...
               y <= halfWindowSize || y >= rows-halfWindowSize
                continue;
            end
            
            % Extract window in img2
            window2 = img2(y-halfWindowSize:y+halfWindowSize, ...
                           x-halfWindowSize:x+halfWindowSize, :);
            
            % Extract window in img1 centered at pts1
            window1 = img1(pts1(i,2)-halfWindowSize:pts1(i,2)+halfWindowSize, ...
                           pts1(i,1)-halfWindowSize:pts1(i,1)+halfWindowSize, :);
            
            % Compute SSD
            SSD = sum((window1(:) - window2(:)).^2);
            
            % Update minimum SSD and its index
            if SSD < minSSD
                minSSD = SSD;
                bestIndex = j;
            end
        end
        
        % Update pts2 with the best matching point
        pts2(i,:) = [X2(bestIndex) Y2(bestIndex)];
    end
end



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


