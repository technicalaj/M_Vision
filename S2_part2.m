

% New extrinsic parameters (assuming the world origin doesn't change)
R1p = R_new;
R2p = R_new;
t1p = -R1p * opticalCenterNew;
t2p = -R2p * opticalCenterNew;
K1p = K_new;
K2p = K_new;
M1 = M1; % Rectification matrix for the first image
M2 = M2; % Rectification matrix for the second image

dispMap = get_disparity(img1_rect, img2_rect, max_disp, win_size);
figure; imshow(dispMap, []); colormap('jet'); colorbar;
title('Disparity Map');

function dispMap = get_disparity(img1_rect, img2_rect, max_disp, win_size)
    [rows, cols] = size(img1_rect);
    dispMap = zeros(rows, cols);
    padSize = floor(win_size / 2);
    
    % Pad images to avoid border issues
    im1_padded = padarray(img1_rect, [padSize padSize], 0, 'both');
    im2_padded = padarray(img2_rect, [padSize padSize], 0, 'both');
    
    for y = 1 + padSize:rows + padSize
        for x = 1 + padSize:cols + padSize
            min_ssd = inf;
            best_disp = 0;
            
            for d = 0:min(max_disp, cols + padSize - x)
                dx = max(1, x-padSize-d);  % Ensure index is not less than 1
                dx_max = min(cols + 2*padSize, x+padSize-d);  % Ensure index does not exceed the image width including padding

                % Extract blocks from both images, ensuring they have matching dimensions
                block1 = im1_padded(y-padSize:y+padSize, x-padSize:x+padSize);
                block2 = im2_padded(y-padSize:y+padSize, dx:dx_max);
                
                % Adjust block1 to match block2's size for fair SSD calculation
                if size(block1, 2) > size(block2, 2)
                    block1 = block1(:, 1:size(block2, 2));
                elseif size(block2, 2) > size(block1, 2)
                    block2 = block2(:, 1:size(block1, 2));
                end
                
                % Compute SSD
                ssd = sum(sum((block1 - block2).^2));
                
                % Update best disparity if current SSD is lower
                if ssd < min_ssd
                    min_ssd = ssd;
                    best_disp = d;
                end
            end
            
            % Update disparity map
            dispMap(y-padSize, x-padSize) = best_disp;
        end
    end
end


