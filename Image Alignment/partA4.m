%% PartA4
clear all; clc; close all;

v1H = VideoReader('video1_high.avi');
v1L = VideoReader('video1_low.avi');
v2H = VideoReader('video2_high.avi');
v2L = VideoReader('video2_low.avi');

vs = input('Select video: (1) video1_high, (2) video1_low, (3) video2_high, (4) video2_low : ');
iter = 100;

if (vs == 1)
    for templ_num = 10:10:100
        [results results_lk MSE,rho,MSELK] = ecc_lk_alignment(read(v1H, templ_num), read(v1H, 1), 1, iter, 'affine', eye(2, 3));
        
        figure(4)        
        subplot(5, 2, templ_num/10);
        hold on
        stem(20*log10(255./MSELK),'r');
        stem(20*log10(255./MSE),'k');
        str = sprintf('PSNR ERRORs Image %d', templ_num); title(str);
        hold off
    end
    
elseif (vs == 2)
    for templ_num = 10:10:100
        [results results_lk MSE,rho,MSELK] = ecc_lk_alignment(read(v1L, templ_num), read(v1L, 1), 1, iter, 'affine', eye(2, 3));
        
        figure(4)        
        subplot(5, 2, templ_num/10);
        hold on
        stem(20*log10(255./MSELK),'r');
        stem(20*log10(255./MSE),'k');
        str = sprintf('PSNR ERRORs Image %d', templ_num); title(str);
        hold off    
    end
        
elseif (vs == 3)
    for templ_num = 10:10:100
        [results results_lk MSE,rho,MSELK] = ecc_lk_alignment(read(v2H, templ_num), read(v2H, 1), 1, iter, 'affine', eye(2, 3));
        
        figure(4)        
        subplot(5, 2, templ_num/10);
        hold on
        stem(20*log10(255./MSELK),'r');
        stem(20*log10(255./MSE),'k');
        str = sprintf('PSNR ERRORs Image %d', templ_num); title(str);
        hold off
    end
    
elseif (vs == 4)
    for templ_num = 10:10:100
        [results results_lk MSE,rho,MSELK] = ecc_lk_alignment(read(v2L, templ_num), read(v2L, 1), 1, iter, 'affine', eye(2, 3));
        
        figure(4)        
        subplot(5, 2, templ_num/10);
        hold on
        stem(20*log10(255./MSELK),'r');
        stem(20*log10(255./MSE),'k');
        str = sprintf('PSNR ERRORs Image %d', templ_num); title(str);
        hold off
    end
    
else
    disp('Choose a number between 1-4 !');
end
legend('LK', 'ECC');