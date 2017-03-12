%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Shared pudding       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

BG(1:256, 1:360, 1:3, 1:17) = 0; % Init the 17 256x360x3 templates 

A = im2double(imread('pudding', 'png')); % Read image

cnt = 1; % Init counter
sz_old = 308; 
step = 0; % Init step

for a = -0.2 : 0.025 : 0.2
    tf = affine2d([1 0 0;
                   a 1 0;
                   0 0 1]); % Generate transformation
    
    AT = imwarp(A,tf); % Apply transformation
    sz_new = size(AT); % Get size of transformed image
    
    step = step + (sz_old - sz_new(1, 2)); % Update step
    sz_old = sz_new(1, 2); 
    
    if(a <= 0) % Cases a = -0.2 to 0
       BG(1:256, 53:(360-step), 1:3, cnt) = AT;
    else % Cases a = 0.025 to 0.2
      BG(1:256, (1 + step):308, 1:3, cnt) = AT;
    end
    BG(BG == 0) = 255;
    cnt = cnt + 1; % Increase counter
end

% Create multiframe indexed image C of matrices BG
C = cat(4, BG(:, :, :, end:-1:1), BG(:, :, :, 2:end),...
           BG(:, :, :, (end-1):-1:1), BG(:, :, :, 2:end));


mov = immovie(C); % Create movie structure array
save('shared_pudding.mat','mov');
implay(mov); % Use Video Viewer to show the movie structure array