% Example script for Part 4: Implement Triangulation

%% Preliminary: Load or Define All Necessary Data


%% Decompose the Essential Matrix to Obtain R and t
[U, ~, V] = svd(E);
W = [0 -1 0; 1 0 0; 0 0 1];
R_options = {U*W*V', U*W'*V'};
t_options = {U(:,3), -U(:,3)};

% Assuming the correct R and t are the first options
% This is a simplification, in practice, you should verify which combination
% works correctly with your scene.
R = R_options{1};
t = t_options{1};

%% Construct Camera Projection Matrices P1 and P2
P1 = K1 * [eye(3), zeros(3, 1)];
P2 = K2 * [R, t];

%% Triangulate Points
% Triangulate points using the Direct Linear Transform (DLT) method
pts3D = triangulatePoints(P1, P2, pts1, pts2);

%% Calculate Reprojection Error
reprojectionError = calculateReprojectionError(P1, P2, pts1, pts2, pts3D);

%% Display Results
disp('Triangulated 3D Points:');
disp(pts3D);

disp('Average Reprojection Error:');
disp(reprojectionError);

%% Helper Functions



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
