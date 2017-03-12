%% PartB 2
clear all, clc;

load('img.mat');

imgv = zeros(size(img, 1), size(img, 2), 10);
for i = 1:10
   imgv(:, : , i) = img; 
end

open('empty.mdl');