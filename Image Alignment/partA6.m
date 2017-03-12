%% PartA6
clear all; clc; close all;

v1H = VideoReader('video1_high.avi');

mean_sq_er_ecc = 0;
mean_sq_er_lk = 0;

iter = 50; % Function iterations
var = 4; % Gaussian grey levels
a = 18^(1/3); % Uniform gray levels

for i = 1:100
% Insert Gaussian noise
noise_image = imnoise(read(v1H, 10),'gaussian', 0, var);

% Insert Uniform noise
% noise_image = double(read(v1H, 10)) + 2*(rand(256,256)-.5)*a;

[results results_lk MSE,rho,MSELK] = ecc_lk_alignment(noise_image, read(v1H, 1), 1, iter, 'affine', eye(2, 3));

mean_sq_er_ecc = mean_sq_er_ecc + MSE(iter);
mean_sq_er_lk = mean_sq_er_lk + MSELK(iter);
end

% Mean square errors for the 100 executions
mean_sq_er_ecc = mean_sq_er_ecc / i; 
mean_sq_er_lk = mean_sq_er_lk / i;