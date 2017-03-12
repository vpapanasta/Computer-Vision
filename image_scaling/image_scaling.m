%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Image scaling        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

% Define 2018x1600x3 zero table
synthesis(1:2018, 1:1600, 1:3) = 0;

MW = im2double(imread('milkyway', 'jpg')); % 538x800x3
sz = size(MW);
image(MW);
title('Milkyway');

% x2 scaling transformation
tform_1 = affine2d([2 0 0;
                    0 2 0;
                    0 0 1]);
% /2 scaling transformation
tform_2 = affine2d([0.5 0 0;
                    0 0.5 0;
                    0 0 1]);
% /4 scaling transformation
tform_3 = affine2d([0.25 0 0;
                    0 0.25 0;
                    0 0 1]);

MW_1 = im2double(imwarp(MW,tform_1)); % 1076x1600x3
sz1 = size(MW_1);

MW_2 = im2double(imwarp(MW,tform_2)); % 269x400x3
sz2 = size(MW_2);

MW_3 = im2double(imwarp(MW,tform_3)); % 135x200x3
sz3 = size(MW_3);

% Create synthesis of scaled images
synthesis(1:sz1(1, 1), 1:sz1(1, 2), :) = MW_1;

synthesis((sz1(1, 1) + 1):(sz1(1, 1) + sz(1, 1)), 1:sz(1, 2), :) = MW;
synthesis((sz1(1, 1) + 1):(sz1(1, 1) + sz(1, 1)), (sz(1, 2) + 1):sz1(1, 2) , :) = MW;

synthesis((sz1(1, 1) + sz(1, 1) + 1):(sz1(1, 1) + sz(1, 1) + sz2(1, 1)), 1:sz2(1, 2), :) = MW_2;
for i = 1:3
    synthesis((sz1(1, 1) + sz(1, 1) + 1):(sz1(1, 1) + sz(1, 1) + sz2(1, 1)), (i*sz2(1, 2) + 1):((i+1)*sz2(1, 2)), :) = MW_2;
end

synthesis((sz1(1, 1) + sz(1, 1) + sz2(1, 1) + 1):(sz1(1, 1) + sz(1, 1) + sz2(1, 1) + sz3(1, 1)), 1:sz3(1, 2), :) = MW_3;
for i = 1:7
    synthesis((sz1(1, 1) + sz(1, 1) + sz2(1, 1) + 1):(sz1(1, 1) + sz(1, 1) + sz2(1, 1) + sz3(1, 1)), (i*sz3(1, 2) + 1):((i+1)*sz3(1, 2)), :) = MW_3;
end

% Plot synthesis of scaled images
figure;
image(synthesis);
title('Milkyway Synthesis');