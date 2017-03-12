function [results results_lk MSE,rho,MSELK]= ecc_lk(image, template, levels, noi, transform, delta_p_init)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ECC image alignment algorithm
%RESULTS = ECC(IMAGE, TEMPLATE, LEVELS, NOI, TRANSFORM, DELTA_P_INIT)
%
% This m-file implements the ECC image alignment algorith which is
% presented in the paper "G.D.Evangelidis, E.Z.Psarakis, Parametric Image Alignment
% using Enhanced Correlation Coefficient.IEEE Trans. on PAMI, vol.30, no.10, 2008"
% The code outline follows at some extent the outline of Matlab code released
% by Baker-et-al for "Lucas-Kanade 20 Years On" project of CMU.
% ------------------
% Input variables:
% IMAGE:        the profile needs to be warped in order to be similar to TEMPLATE,
% TEMPLATE:     the profile needs to be reached,
% NOI:          the number of iterations per level; the algorithm is executed
%               (NOI-1) times
% LEVELS:       the number of levels in pyramid scheme (set LEVELS=1 for a
%               non pyramid implementation), the level index 1
%               corresponds to the level with the highest image resolution
% TRANSFORM:    the type of adopted transform, accepted strings: 'affine','homography'
% DELTA_P_INIT: the initial transformation matrix for original images (optional); The identity
%               transformation is the default value (see 'trasnform initialization'
%               subroutine). In case of affine transform, DELTA_P_INIT must be a
%               2x3 matrix, while in homography case, it must be a 3x3 matrix.
%
% Output:
% RESULTS:   A struct of size LEVELSxNOI with the following fileds:
%
% RESULTS().warp:              the warp needs to be applied in image at each level-iteration,
% RESULTS().rho:               the correlation coefficient value at each level-iteration,
% RESULTS(LEVELS,NOI).image:   the final warped image which is similar to TEMPLATE.
%
% The first stored .warp and .rho values are due to the initialization. In
% case of pour final alignment results check the initialization of the
% algorithm and/or overlap of the images.
% -------------------
% $ Ver: 1.0.0, 1/3/2010,  released by Georgios D. Evangelidis, Fraunhofer IAIS.
% For any comment, please contact georgios.evangelidis@iais.fraunhofer.de
% or evagelid@ceid.upatras.gr
%
% This software is for research purposes only. In any case, please cite the above paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all
tic;

% If plot_flag=1, initial and final images are plotted at the end of
% execution
plot_flag=1;

break_flag=0;
if nargin<5
    error('-> Not enough input arguments');
end

if ~(strcmp(transform,'affine')||strcmp(transform,'homography'))
    error('-> Not a valid transform string')
end

template = double(template);
template=255*(template-min(template(:)))/max(template(:)-min(template(:)));
image = double(image);
image=255*(image-min(image(:)))/max(image(:)-min(image(:)));
%% pyramid images
% The following for-loop creates pyramid images with varying names
% im_1,im_2,...im_levels and temp_1,temp_2,...,temp_levels. The
% images im_1, temp_1 are the images with the highest resoltuion

% Smoothing of original images
f = fspecial('gaussian',[5 5],.5);
% temp_1 = filter2(f,template);
% im_1 = filter2(f,image);
temp_1=template;
im_1=image;
for nol=2:levels
    eval(['im_' num2str(nol) '=imresize(im_' num2str(nol-1) ',.5);'])
    eval(['temp_' num2str(nol) '=imresize(temp_' num2str(nol-1) ',.5);'])
end
%LK
im_lk=im_1;

%% transform initialization

% In case of affine transform the mapping function (see the above paper) has the form
%   phi(x;p)= [1+p1, p3, p5;
%              p2, 1+p4, p6]
%  and
%  delta_p = [p1, p3, p5;
%             p2, p4, p6]

% In case of homography transform the mapping function (see the above paper) has the form
%   phi(x;p)= [1+p1, p4, p7;
%              p2, 1+p5, p8;
%              p3,   p6,  1]
%  and
%  delta_p = [p1, p4, p7;
%             p2, p5, p8;
%             p3, p6,  1]


if strcmp(transform,'affine')
    nop=6; %number of parameters
    if nargin==5;
        warp=zeros(3,3);
        warp_lk=zeros(3,3);
    else
        if (size(delta_p_init,1)~=2)|(size(delta_p_init,2)~=3)
            error('-> In affine case the size of initialization matrix must be 2x3');
        else
            warp=[delta_p_init;zeros(1,3)];
            warp_lk=[delta_p_init;zeros(1,3)];
        end
    end
end

if strcmp(transform,'homography')
    nop=8; %number of parameteres
    if nargin==5;
        warp=zeros(3,3);
    else
        if (size(delta_p_init,1)~=3)|(size(delta_p_init,2)~=3)
            error('-> In homography case the size of initialization matrix must be 3x3');
        else
            warp=delta_p_init;
        end
    end
end

% in case of pyramid implementation, the initial transformation must be
% appropriately modified
for ii=1:levels-1
    warp=next_level(warp, transform, 0);
end

%% Run ECC algorithm for each level of pyramid
for nol=levels:-1:1

    eval(['im=im_' num2str(nol) ';'])

    [vx,vy]=gradient(im);
    [vx_lk,vy_lk]=gradient(im_lk);
    eval(['temp=temp_' num2str(nol) ';'])

    [A,B]=size(temp);
    % Warning for tiny images
    if prod([A,B])<400
        disp('-> ECC Warning: The size of images in high levels is quite small and the results may be affected.');
        disp('-> To avoid this case try fewer levels or bigger images.');
        disp('-> Press any button to continue.')
        pause
    end


    %Define the rectangular Region of Interest by nx and ny (you can modify the ROI).
    %Here we just ignore image margins. Margin is equal to 5 percent of the mean of [height,width].

    m0=mean([A,B]);
    margin=floor(m0*.05/(2^(nol-1)));

    nx=margin+1:B-margin;
    ny=margin+1:A-margin;
    temp=double(temp(ny,nx,:));

    temp=temp-mean(temp(:)); % zero-mean image; is useful for brithness change compensation, otherwise you can comment this line
    n_temp=norm(temp(:));
    %temp=temp/n_temp;     % this normalization does not affect the results. The closed-form solution below is invariant to this normalization.

%LK
im_lk=im;
    %% ECC, Forwards Additive Algorithm -------------------------------
    for i=1:noi

        disp(['Level: ' num2str(nol) ', Iteration: ' num2str(i)])
        %Image interpolation method
        str='bilinear'; % bilinear interpolation
        wim = spatial_interp(im, warp, str, transform, nx, ny);
        wim = wim-mean(wim(:));% zero-mean image; is useful for brithness change compensation, otherwise you can comment this line
        wim_lk = spatial_interp(im_lk, warp_lk, str, transform, nx, ny);
        wim_lk = wim_lk-mean(wim_lk(:));% zero-mean image; is useful for brithness change compensation, otherwise you can comment this line
        %Save current transform
        if strcmp(transform,'affine')
            results(nol,i).warp = warp(1:2,:);
            results_lk(nol,i).warp = warp_lk(1:2,:);
        else
            results(nol,i).warp = warp;
        end

        results(nol,i).rho = dot(temp(:),wim(:)) / n_temp / norm(wim(:));
        results_lk(nol,i).rho = dot(temp(:),wim_lk(:)) / n_temp / norm(wim_lk(:));

        image_error=temp(:)-wim(:);
        image_error_lk=temp(:)-wim_lk(:);

        MSE(i)=sqrt(mean2(image_error.^2));
        rho(i)=dot(temp(:),wim(:)) / n_temp / norm(wim(:));
        MSELK(i)=sqrt(mean2(image_error_lk.^2));
        if (i == noi) % the algorithm is executed (noi-1) times
            break;
        end

        % Gradient Image interpolation (warped gradients)
        wvx = spatial_interp(vx, warp, str, transform, nx, ny);
        wvy = spatial_interp(vy, warp, str, transform, nx, ny);


        % Compute the jacobian of warp transform
        J = warp_jacobian(nx, ny, warp, transform);

        % Compute the jacobian of wim wrt parameters (matrix G in paper)
        G = image_jacobian(wvx, wvy, J, nop);

        % Compute Hessian and its inverse
        C= G' * G;
        
        con=cond(C);

        if con>1.0e+15
            disp('->ECC Warning: Badly conditioned Hessian matrix. Check the initialization or the overlap of images.')
        end    
        i_C = inv(C);
        
        % Compute projections of image vectors into G-space
        Gt = G' * temp(:);
        Gw = G' * wim(:);

%LK 
% Gradient Image interpolation (warped gradients)
        wvx_lk = spatial_interp(vx_lk, warp_lk, str, transform, nx, ny);
        wvy_lk = spatial_interp(vy_lk, warp_lk, str, transform, nx, ny);


        % Compute the jacobian of warp transform
        J_lk = warp_jacobian(nx, ny, warp_lk, transform);

        % Compute the jacobian of wim wrt parameters (matrix G in paper)
        G_lk = image_jacobian(wvx_lk, wvy_lk, J_lk, nop);

        % Compute Hessian and its inverse
        C_lk= G_lk' * G_lk;
        
        con=cond(C_lk);

        if con>1.0e+15
            disp('->LK Warning: Badly conditioned Hessian matrix. Check the initialization or the overlap of images.')
        end    
        
        
        %% ECC closed form solution

        % Compute lambda parameter
        num = (norm(wim(:))^2 - Gw' * i_C * Gw);
        den = (dot(temp(:),wim(:)) - Gt' * i_C * Gw);
        lambda = num / den;
        T_proj=Gt'*i_C*Gt/(norm(temp(:))*norm(temp(:)));
        W_proj=Gw'*i_C*Gw/(norm(wim(:))*norm(wim(:)));
        TW_proj=Gw'*i_C*Gt/(norm(temp(:))*norm(wim(:)));
        E_proj=(Gw/norm(wim(:))-Gt/norm(temp(:)))'*i_C*(Gw/norm(wim(:))-Gt/norm(temp(:)));
        
        % Compute error vector
        imerror = lambda * temp -wim;
        
        % Compute the projection of error vector into Jacobian G
        Ge = G' * imerror(:);

        % Compute the optimum parameter correction vector
        delta_p = i_C * Ge;


        if (sum(isnan(delta_p)))>0 %Hessian is close to singular
            disp(['-> Algorithms stopped at ' num2str(i) '-th iteration of ' num2str(nol) '-th level due to bad condition of Hessian matrix.']);
            disp(['-> Final results are stored at results(' num2str(nol) ',' num2str(i) ').warp and results(' num2str(nol) ',' num2str(i) ').image. Note that ']);
            disp('-> final warp and image have been modified with respect to original size of images.');

            break_flag=1;
            break;
        end

        % Update parmaters
        warp = param_update(warp, delta_p, transform);


%% LK closed form solution
% Compute parameters
        a_lk =temp(:)'/n_temp*wim_lk(:);
        b_lk=G_lk'*temp(:)/n_temp;
        
        % Compute error vector
        imerror_lk =a_lk*temp(:)/n_temp-wim_lk(:);

        % Compute the projection of error vector into Jacobian G
        Ge_lk = G_lk' * imerror_lk(:);
        i_C_lk=inv(G_lk'*G_lk-b_lk*b_lk');
        % Compute the optimum parameter correction vector
        delta_p_lk = i_C_lk * Ge_lk;


        if (sum(isnan(delta_p_lk)))>0 %Hessian is close to singular
            disp(['-> Algorithms stopped at ' num2str(i) '-th iteration of ' num2str(nol) '-th level due to bad condition of Hessian matrix.']);
            disp(['-> Final results are stored at results(' num2str(nol) ',' num2str(i) ').warp and results(' num2str(nol) ',' num2str(i) ').image. Note that ']);
            disp('-> final warp and image have been modified with respect to original size of images.');

            break_flag=1;
            break;
        end

        % Update parmaters
        warp_lk = param_update(warp_lk, delta_p_lk, transform);


    end


    if break_flag==1
        break;
    end

    % modify the parameteres appropriately for next pyramid level
    if (nol>1)&(break_flag==0)
        warp = next_level(warp, transform,1);
        warp_lk = next_level(warp_lk, transform,1);
    end

end

toc

if break_flag==1 % this conditional part is only executed when algorithm stops due to Hessian singularity
    for jj=1:nol-1
        warp = next_level(warp, transform,1);
        m0=2*m0;
    end
    margin=floor(m0*.05);
    nx=margin+1:size(template,2)-margin;
    ny=margin+1:size(template,1)-margin;
    results(nol,i).warp=warp;
    results_lk(nol,i).warp=warp_lk;
end

% store the final warped image
results(nol,i).image = spatial_interp(image, results(nol,i).warp, str, transform, nx, ny);
results_lk(nol,i).image = spatial_interp(image, results_lk(nol,i).warp, str, transform, nx, ny);
% project ROI corners through final warp
ROI_corners=[nx(1) nx(1) nx(end) nx(end);...
    ny(1) ny(end) ny(1) ny(end)];

Mat=results(nol,i).warp;
Mat(1:2,1:2)=Mat(1:2,1:2)+eye(2);
wROI_corners=Mat*[ROI_corners;ones(1,4)];

Mat_lk=results_lk(nol,i).warp;
Mat_lk(1:2,1:2)=Mat_lk(1:2,1:2)+eye(2);
wROI_corners_lk=Mat_lk*[ROI_corners;ones(1,4)];

if strcmp(transform,'homography')
    wROI_corners=wROI_corners./repmat(wROI_corners(3,:),3,1);
end


if plot_flag==1
    figure(1)
    % plot images for highest-resolution level of pyramid
    subplot(2,2,1)
    imshow(uint8(template))
    hold on
    line([nx(1) nx(end)],[ny(1) ny(1)],'Color','m')
    line([nx(end) nx(end)],[ny(1) ny(end)],'Color','m')
    line([nx(1) nx(end)],[ny(end) ny(end)],'Color','m')
    line([nx(1) nx(1)],[ny(1) ny(end)],'Color','m')
    hold off
    title('Template with marked ROI')
    axis on

    subplot(2,2,2)
    imshow(uint8(image))
    hold on
    line([wROI_corners(1,1) wROI_corners(1,3)],[wROI_corners(2,1) wROI_corners(2,3)],'Color','m')
    line([wROI_corners(1,3) wROI_corners(1,4)],[wROI_corners(2,3) wROI_corners(2,4)],'Color','m')
    line([wROI_corners(1,2) wROI_corners(1,4)],[wROI_corners(2,2) wROI_corners(2,4)],'Color','m')
    line([wROI_corners(1,1) wROI_corners(1,2)],[wROI_corners(2,1) wROI_corners(2,2)],'Color','m')
    hold off
    title('Input image with warped ROI')
    axis on

    subplot(2,2,3)
    imshow(uint8(results(nol,i).image))
    title('Warped image')
    axis on


    image_error=(results(nol,i).image-template(ny,nx));
    image_error=image_error-min(image_error(:));
    image_error=image_error./max(image_error(:))*255;
    subplot(2,2,4)
    imshow(uint8(image_error))
    title('Error image')
    axis on
    %LK
    figure(2)
    subplot(2,2,1)
    imshow(uint8(template))
    hold on
    line([nx(1) nx(end)],[ny(1) ny(1)],'Color','m')
    line([nx(end) nx(end)],[ny(1) ny(end)],'Color','m')
    line([nx(1) nx(end)],[ny(end) ny(end)],'Color','m')
    line([nx(1) nx(1)],[ny(1) ny(end)],'Color','m')
    hold off
    title('Template with marked ROI')
    axis on
    
    subplot(2,2,2)
    imshow(uint8(image))
    hold on
    line([wROI_corners_lk(1,1) wROI_corners_lk(1,3)],[wROI_corners_lk(2,1) wROI_corners_lk(2,3)],'Color','m')
    line([wROI_corners_lk(1,3) wROI_corners_lk(1,4)],[wROI_corners_lk(2,3) wROI_corners_lk(2,4)],'Color','m')
    line([wROI_corners_lk(1,2) wROI_corners_lk(1,4)],[wROI_corners_lk(2,2) wROI_corners_lk(2,4)],'Color','m')
    line([wROI_corners_lk(1,1) wROI_corners_lk(1,2)],[wROI_corners_lk(2,1) wROI_corners_lk(2,2)],'Color','m')
    hold off
    title('Input image with warped ROI')
    axis on

    subplot(2,2,3)
    imshow(uint8(results_lk(nol,i).image))
    title('Warped image')
    axis on


    image_error=(results_lk(nol,i).image-template(ny,nx));
    image_error=image_error-min(image_error(:));
    image_error=image_error./max(image_error(:))*255;
    subplot(2,2,4)
    imshow(uint8(image_error))
    title('Error image')
    axis on
%ERRORS
figure(3)
    title('PSNR ERRORs:(Red:LK)-(Black:ECC)')
    hold on
stem(20*log10(255./MSELK),'r');
stem(20*log10(255./MSE),'k');
hold off
end

