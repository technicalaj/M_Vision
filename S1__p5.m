%% Step 1: Initialization and Load Data
% Ensure these files are in your MATLAB path or current directory
load('data/intrinsics.mat'); % Contains K1, K2
load('data/some_corresp.mat'); % Contains pts1, pts2
load('data/temple_coords.mat'); % Make sure this contains points from img1 for 3D reconstruction


% Define image paths for the stereo images
img1Path = 'data/im1.png'; % Replace with actual path
img2Path = 'data/im2.png'; % Replace with actual path
img1 = imread('data/im1.png');
img2 = imread('data/im2.png');

% ---> You might need to rename 'points' based on the actual variable name in 'temple_coords.mat'
% If 'points' variable does not exist, replace 'points' with the correct variable name from the .mat file
pointsFromImg1 = points; % Change 'points' to the actual variable name if different


%% Step 2: Compute or Load the Fundamental Matrix (F)
% Assuming 'F' has already been computed in previous steps and is available in the workspace.
% If not, recompute it using your eight_point function or load it from where it's saved.

%% Step 3: Find Epipolar Correspondences
% Adapt the 'findEpipolarCorrespondences' function call as per your function definition
% Ensure this function is defined in your path
pts2 = epipolar_correspondences(img1, img2, F, pts1);


%% Step 4: Compute the Essential Matrix (E) from F
E = K2' * F * K1; % Make sure K1 and K2 are correctly loaded from 'intrinsics.mat'

%% Step 5: Decompose E to Obtain the Correct [R|t] Configuration
[U, ~, V] = svd(E);
    W = [0 -1 0; 1 0 0; 0 0 1]; % Auxiliary matrix for SVD decomposition
    
    % Possible rotation and translation configurations
    R1 = U * W * V';
    R2 = U * W' * V';
    if det(R1) < 0, R1 = -R1; end
    if det(R2) < 0, R2 = -R2; end
    T = U(:, 3);

    configs = [struct('R', R1, 't', T), struct('R', R2, 't', T), ...
               struct('R', R1, 't', -T), struct('R', R2, 't', -T)];

    % Initialize variables to keep track of the best configuration
    maxInFront = 0;
    R_final = [];
    t_final = [];

    for i = 1:length(configs)
        R = configs(i).R;
        t = configs(i).t;

        % Construct camera matrices for this configuration
        P1 = K1 * [eye(3), zeros(3, 1)];
        P2 = K2 * [R, t];

        % Triangulate points
        pts3D = triangulatePoints(P1, P2, pts1, pts2);

        % Check how many points are in front of both cameras
        inFront1 = sum(pts3D(:, 3) > 0); % For camera 1, Z should be positive
        transformedPts = (R * pts3D' + t)'; % Transform points to camera 2's coordinate frame
        inFront2 = sum(transformedPts(:, 3) > 0); % For camera 2, Z should also be positive

        if inFront1 + inFront2 > maxInFront
            maxInFront = inFront1 + inFront2;
            R_final = R;
            t_final = t;
        end
    end




%% Step 8: Visualization of 3D Reconstruction
figure;
scatter3(pts3D(:,1), pts3D(:,2), pts3D(:,3), 'filled');
title('3D Reconstruction of the Temple');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on; rotate3d on;


% Assuming R1, R2, t1, t2 are your final extrinsic parameters from Section 1 Part 5
save('data/extrinsics.mat', 'R1', 'R2', 't1', 't2');


%% Helper Functions Definitions
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




function pts3D = triangulatePoints(P1, P2, pts1, pts2)

    pts1 = double(pts1); % Convert pts1 to double
    pts2 = double(pts2); % Convert pts2 to double
    N = size(pts1, 1);
    pts3D = zeros(N, 3);
    for i = 1:N
        A = [pts1(i,2)*P1(3,:) - P1(2,:);
             P1(1,:) - pts1(i,1)*P1(3,:);
             pts2(i,2)*P2(3,:) - P2(2,:);
             P2(1,:) - pts2(i,1)*P2(3,:)];
        [~, ~, V] = svd(A);
        X = V(:,end);
        pts3D(i,:) = X(1:3)'/X(4); % Convert X to homogeneous coordinates and normalize
    end
end

function error = calculateReprojectionError(P1, P2, pts1, pts2, pts3D)
    % Ensure all inputs are double
    pts1 = double(pts1);
    pts2 = double(pts2);
    pts3D_hom = [double(pts3D), ones(size(pts3D, 1), 1)]'; % Ensure pts3D is double and make homogeneous

    % Project 3D points back to 2D using the camera matrices
    pts1_proj_hom = P1 * pts3D_hom;
    pts2_proj_hom = P2 * pts3D_hom;

    % Convert from homogeneous to Cartesian coordinates
    pts1_proj = pts1_proj_hom(1:2, :) ./ pts1_proj_hom(3, :);
    pts2_proj = pts2_proj_hom(1:2, :) ./ pts2_proj_hom(3, :);

    % Calculate the squared differences
    diffs1 = (pts1' - pts1_proj).^2;
    diffs2 = (pts2' - pts2_proj).^2;

    % Compute the root sum squared error for each point and then average
    error1 = sqrt(sum(diffs1, 1));
    error2 = sqrt(sum(diffs2, 1));
    error = mean([error1, error2]);
end









