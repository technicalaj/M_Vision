% Assuming dispMap, K1, and other variables from Part 2 are already in the workspace

% Baseline distance between the camera centers
% This needs to be adjusted based on your stereo camera setup
baseline = 0.1;  % Example value in meters or the units used in your camera matrix

% Calculate the depth map from the disparity map
depthMap = get_depth(dispMap, K1, baseline);





% Cap the maximum depth for better visualization if needed
maxDepthForVisualization = 152.04; % or any value you find suitable based on your application
depthMapVis = min(depthMap, maxDepthForVisualization);

% Normalize depth map for better visualization
depthMapNorm = (depthMapVis - min(depthMapVis(:))) / (maxDepthForVisualization - min(depthMapVis(:)));

figure;
imshow(depthMapNorm, []);
colormap('jet');
colorbar;
title('Depth Map');

function depthMap = get_depth(dispMap, K1, baseline)   

% The focal length in pixels (assuming square pixels)
    f = K1(1, 1);  % Assuming fx = fy for simplicity

    % Avoid division by zero for pixels with zero disparity
    dispMap(dispMap == 0) = inf;  % Set zero disparities to infinity
    
    % Convert disparity to depth
    depthMap = (f * baseline) ./ dispMap;
    
    % Set infinite depth values (originally zero disparity) to NaN or another identifier
    depthMap(isinf(depthMap)) = NaN;  % Mark as NaN or use a specific value to indicate no depth
end