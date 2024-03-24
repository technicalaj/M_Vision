% Load your data and images
load('data/intrinsics.mat'); % Contains K1, K2
img1 = imread('data/im1.png');
img2 = imread('data/im2.png');



% Compute the camera centers from rotation and translation matrices
C1 = -R1'*t1; % Camera 1 center
C2 = -R2'*t2; % Camera 2 center

% Assuming both cameras have the same intrinsic parameters (for simplicity)
% and your actual setup might require using K1 and K2 separately
K_new = K1; 

% Example transformation matrices for rectification, replace with actual calculations
% Here you need to calculate M1 and M2 based on your rectification process
M1 = eye(3); % Placeholder for the actual transformation matrix for img1
M2 = eye(3); % Placeholder for the actual transformation matrix for img2

% Transform corners of the original images to find the bounding box in rectified space
[h1, w1, ~] = size(img1);
[h2, w2, ~] = size(img2);
corners1 = [1, 1; w1, 1; 1, h1; w1, h1];
corners2 = [1, 1; w2, 1; 1, h2; w2, h2];

% Apply the rectification transformations to the corners
rectifiedCorners1 = transformPointsForward(projective2d(M1'), corners1);
rectifiedCorners2 = transformPointsForward(projective2d(M2'), corners2);

% Determine bounding box in rectified space
allCorners = [rectifiedCorners1; rectifiedCorners2]; % Combine for both images
xLimits = [min(allCorners(:,1)), max(allCorners(:,1))];
yLimits = [min(allCorners(:,2)), max(allCorners(:,2))];

% Expand the limits to ensure no image content is cut off, adjust as needed
padding = 50; % Example padding to ensure no content is lost
xLimits = xLimits + [-padding, padding];
yLimits = yLimits + [-padding, padding];

% Set output view using calculated limits
outputView = imref2d([round(diff(yLimits)), round(diff(xLimits))], ...
                     [min(xLimits), max(xLimits)], [min(yLimits), max(yLimits)]);

% Rectify the images using the calculated output view
img1_rect = imwarp(img1, projective2d(M1'), 'OutputView', outputView);
img2_rect = imwarp(img2, projective2d(M2'), 'OutputView', outputView);

% Display the rectified images for verification
figure;
subplot(1,2,1); imshow(img1_rect); title('Rectified Image 1');
subplot(1,2,2); imshow(img2_rect); title('Rectified Image 2');
