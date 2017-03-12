%%%%%%%%%%%%%%%%%%%%
%      Ball 2      %
%%%%%%%%%%%%%%%%%%%%
clear all; close all;

interp = 'cubic';

b = imread('ball', 'jpg'); % 
bm = imread('ball_mask', 'jpg'); % 
beach = imread('beach', 'jpg'); % 

backg_size = size(beach); % Size of background image
width = backg_size(1, 1); % Width of background image
height = backg_size(1, 2); % Height of background image 

temp = uint8(zeros(width, height, 3, 63)); % Define the template frames 

b = b - bm; % Sub ball image with the mask
bm = imcomplement(bm); % Reverse mask colors

refresh = beach; % Backup of background matrix

%%%%%%%%% Ball Scaling %%%%%%%%%
scale = affine2d([0.25 0 0;
                  0 0.25 0;
                  0 0 1]); % Generate transformation
% Apply transformation 
[b, b_ref1] = imwarp(b,scale);
[bm, bm_ref1] = imwarp(bm,scale); 

ball_size = size(bm);
wb = ball_size(1, 1);
hb = ball_size(1, 2);

theta = 10; % Define theta
% Init scalars of y-axis translation function
a = 0.2; d = 5; c = 300; 
k = 1; % Init counter

% For each t
for t = 0:0.1:2*pi
    
    %%%%%%%%% Center Translation %%%%%%%%%
    center_trans = affine2d([1 0 0;
                             0 1 0;
                            -hb/2 -wb/2 1]); % Generate transformation
     % Apply transformation          
    [b_t, b_ref] = imwarp(b, b_ref1, center_trans); 
    [bm_t, bm_ref] = imwarp(bm, bm_ref1, center_trans);
  
    %%%%%%%%% Rotation %%%%%%%%%
    rot = affine2d([cosd(theta) -sind(theta) 0;
                    sind(theta) cosd(theta) 0; 0 0 1]); 
    % Apply transformation 
    [b_t, b_ref]  = imwarp(b_t, b_ref, rot, interp, 'FillValues', 0);              
    [bm_t, bm_ref] = imwarp(bm_t, bm_ref, rot, interp, 'FillValues', 0); 
     
    %%%%%%%%% InvCenter Translation %%%%%%%%%
    center_trans = affine2d([1 0 0;
                             0 1 0;
                             hb/2+300 wb/2+300 1]); % Generate transformation
    % Apply transformation          
    [b_t, b_ref] = imwarp(b_t, b_ref, center_trans);
    [bm_t, bm_ref] = imwarp(bm_t, bm_ref, center_trans);
    
      
    %%%%%%%%% Translation %%%%%%%%%
    % Calculate y-axis translation function
    p = 300 - abs(c*cos(d*t)).*exp(-a*t); 
    
    tran = affine2d([1 0 0;
                     0 1 0;
                     p 1 1]); % Generate transformation
    % Apply transformation              
    [b_t, b_ref] = imwarp(b_t, b_ref, tran, interp, 'FillValues', 0);             
    [bm_t, bm_ref] = imwarp(bm_t, bm_ref, tran, interp, 'FillValues', 0);
    
    %%%%%%%%% Ball Scaling %%%%%%%%%
    scale = affine2d([1/k 0 0;
                      0 1/k 0;
                      0 0 1]); % Generate transformation
    % Apply transformation 
    [b_t, b_ref] = imwarp(b_t, b_ref, scale); 
    [bm_t, bm_ref] = imwarp(bm_t, bm_ref, scale);
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(k == 1)
        xx_start = floor(b_ref.XWorldLimits(1)) + 220;
        xx_end = floor(b_ref.XWorldLimits(2)) + 220;
        yy_start = floor(b_ref.YWorldLimits(1)) + 600;
        yy_end = floor(b_ref.YWorldLimits(2)) + 600;
    else
        xx_start = floor(b_ref.XWorldLimits(1))+ 350 + k/2;
        xx_end = floor(b_ref.XWorldLimits(2)) + 350 + k/2;
        yy_start = floor(b_ref.YWorldLimits(1)) + 600;
        yy_end = floor(b_ref.YWorldLimits(2)) + 600;
    end
    
    bt_size = size(b_t);
    
    if(bt_size(1, 1) == bt_size(1, 2)) 
        beach((xx_start + 1):xx_end, (yy_start + 1):yy_end, :) = beach((xx_start + 1):xx_end, (yy_start + 1):yy_end, :) - bm_t;
        beach((xx_start + 1):xx_end, (yy_start + 1):yy_end, :) = beach((xx_start + 1):xx_end, (yy_start + 1):yy_end, :) + b_t;
    elseif(bt_size(1, 1) < bt_size(1, 2))        
        beach((xx_start + 2):xx_end, (yy_start):yy_end, :) = beach((xx_start + 2):xx_end, (yy_start):yy_end, :) - bm_t;
        beach((xx_start + 2):xx_end, (yy_start):yy_end, :) = beach((xx_start + 2):xx_end, (yy_start):yy_end, :) + b_t;
    else        
        beach((xx_start):xx_end, (yy_start + 2):yy_end, :) = beach((xx_start):xx_end, (yy_start + 2):yy_end, :) - bm_t;
        beach((xx_start):xx_end, (yy_start + 2):yy_end, :) = beach((xx_start):xx_end, (yy_start + 2):yy_end, :) + b_t;
    end
    
    temp(:, :, :, k) = beach; % Save new frame
    theta = theta + 10; % Increase angle
    k = k + 1; % Increase counter
    beach = refresh; % Background refresh
    
end

% Create multiframe indexed image C of matrices BG
C = cat(4, temp(:, :, :, 1:end));

mov = immovie(C); % Create movie structure array
save('ball2.mat','mov');
implay(mov, 10);% Use Video Viewer to show the movie structure array