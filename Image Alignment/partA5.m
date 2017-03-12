%% PartA5
clear all; clc; close all;

v1H = VideoReader('video1_high.avi');
iter = 100; % Function iterations
a = 4; b = 10; 

[results results_lk MSE,rho,MSELK] = ecc_lk_alignment(read(v1H, 50), a*read(v1H, 1) + b, 1, iter, 'affine', eye(2, 3));
corr_coef_ecc = results(iter).rho
corr_coef_lk = results_lk(iter).rho

close all;
[results results_lk MSE,rho,MSELK] = ecc_lk_alignment(a*read(v1H, 50) + b, read(v1H, 1), 1, iter, 'affine', eye(2, 3));
corr_coef_ecc = results(iter).rho
corr_coef_lk = results_lk(iter).rho

close all;
[results results_lk MSE,rho,MSELK] = ecc_lk_alignment(a*read(v1H, 50) + b, a*read(v1H, 1) + b, 1, iter, 'affine', eye(2, 3));
corr_coef_ecc = results(iter).rho
corr_coef_lk = results_lk(iter).rho
