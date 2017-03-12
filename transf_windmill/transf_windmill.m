%%%%%%%%%%%%%%%%%%%%
%     Windmill     %
%%%%%%%%%%%%%%%%%%%%
clear all; close all;

interp = '';

interp = input('Please select intepolation method: ', 's');

A = imread('windmill', 'png', 'BackgroundColor', [1 1 1]); % Read windmill image
B = imread('windmill_mask', 'png','BackgroundColor', [1 1 1]); % Read mask image
backg = imread('windmill_back', 'jpeg'); % Read background image

backg_size = size(backg); % Size of background image

width = backg_size(1, 1); % Width of background image
height = backg_size(1, 2); % Height of background image 

temp = uint8(zeros(width, height, 3, 19)); % Define the template frames  

A = A - B; % Sub windmill image with the mask
B = imcomplement(B); % Reverse mask colors

refresh = backg; % Backup of background matrix

cnt = 1; % Define counter

% For angles from 0 to 90 degrees
for theta = 0:5:90 
    % Rotation transform
    tform = affine2d([cosd(theta) -sind(theta) 0;
                      sind(theta) cosd(theta) 0; 0 0 1]); 
    % Apply transformation 
    A_T = imwarp(A,tform, interp, 'FillValues', 0);              
    B_T = imwarp(B,tform, interp, 'FillValues', 0); 
   
    sz = size(B_T); % Size of transformed matrix
    p1 = (960 - sz(1, 1))/2;
    p2 = (1280 - sz(1, 1))/2;
    
    backg((p1+1):(width-p1),(p2+1):(height-p2), :) = backg((p1+1):(width-p1),(p2+1):(height-p2), :) - B_T;
    backg((p1+1):(width-p1),(p2+1):(height-p2), :) = backg((p1+1):(width-p1),(p2+1):(height-p2), :) + A_T;
    
    temp(:, :, :, cnt) = backg; % Save new frame
    backg = refresh; % Background refresh
    cnt = cnt + 1; % Increase counter
end

% Create multiframe indexed image C of matrices BG
C = cat(4, temp(:, :, :, 1:end), temp(:, :, :, 2:end),...
           temp(:, :, :, 2:end), temp(:, :, :, 2:end));


mov = immovie(C); % Create movie structure array
save('transf_windmill.mat','mov');
implay(mov); % Use Video Viewer to show the movie structure array