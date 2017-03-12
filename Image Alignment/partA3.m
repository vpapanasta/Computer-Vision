%% PartA3
clear all; clc; close all;

v1H = VideoReader('video1_high.avi');
v1L = VideoReader('video1_low.avi');
v2H = VideoReader('video2_high.avi');
v2L = VideoReader('video2_low.avi');

vs = input('Select video: (1) video1_high, (2) video1_low, (3) video2_high, (4) video2_low : ');
templ_num = input('Insert template number: ');
iter = 60;

if (vs == 1)
    [results results_lk MSE,rho,MSELK] = ecc_lk_alignment(read(v1H, templ_num), read(v1H, 1), 1, iter, 'affine', eye(2, 3));
    corr_coef_ecc = results(iter).rho
    corr_coef_lk = results_lk(iter).rho
    
    figure(4);
    subplot(1, 2, 1);  
    plot(1:60, MSE);
    title('MSE ECC'); xlabel('Iterations'); ylabel('Mean Square Error');
    subplot(1, 2, 2); 
    plot(1:60, MSELK);
    title('MSE LK'); xlabel('Iterations'); ylabel('Mean Square Error'); 
elseif (vs == 2)
    [results results_lk MSE,rho,MSELK] = ecc_lk_alignment(read(v1L, templ_num), read(v1L, 1), 1, iter, 'affine', eye(2, 3));
    corr_coef_ecc = results(iter).rho
    corr_coef_lk = results_lk(iter).rho
    
    figure(4);
    subplot(1, 2, 1);  
    plot(1:60, MSE);
    title('MSE ECC'); xlabel('Iterations'); ylabel('Mean Square Error');
    subplot(1, 2, 2); 
    plot(1:60, MSELK);
    title('MSE LK'); xlabel('Iterations'); ylabel('Mean Square Error');
elseif (vs == 3)
    [results results_lk MSE,rho,MSELK] = ecc_lk_alignment(read(v2H, templ_num), read(v2H, 1), 1, iter, 'affine', eye(2, 3));
    corr_coef_ecc = results(iter).rho
    corr_coef_lk = results_lk(iter).rho
    
    figure(4);
    subplot(1, 2, 1);  
    plot(1:60, MSE);
    title('MSE ECC'); xlabel('Iterations'); ylabel('Mean Square Error');
    subplot(1, 2, 2); 
    plot(1:60, MSELK);
    title('MSE LK'); xlabel('Iterations'); ylabel('Mean Square Error');
elseif (vs == 4)
    [results results_lk MSE,rho,MSELK] = ecc_lk_alignment(read(v2L, templ_num), read(v2L, 1), 1, iter, 'affine', eye(2, 3));
    corr_coef_ecc = results(iter).rho
    corr_coef_lk = results_lk(iter).rho
    
    figure(4);
    subplot(1, 2, 1);  
    plot(1:60, MSE);
    title('MSE ECC'); xlabel('Iterations'); ylabel('Mean Square Error');
    subplot(1, 2, 2); 
    plot(1:60, MSELK);
    title('MSE LK'); xlabel('Iterations'); ylabel('Mean Square Error');
else
    disp('Choose a number between 1-4 !');
end