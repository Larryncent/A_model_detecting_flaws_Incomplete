% Modified by Xian-Hong Wang June 2023.
% Based on published code TILT_v1.03 by Kelvin Zhang, Arvind Ganesh, February 2011. 
%
% Reference: Linearized Alternating Direction Method with Adaptive Penalty
%            for Fast Solving Transform Invariant Low-rank Texture.
%            Xiang Ren, Zhouchen Lin. Submitted to IJCV.
%
%            TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc.
%            of ACCV, 2010.

% main_affine.m will show you how to run TILT (implemented by Linearized Alternating 
% Direction Method with Adaptive Penalty and warm start techniques) 
% affine on an image, step by step.

clc;
clear;
close all;
task_name='(Result)AFFINE';

%% load the image
[img_name, img_path]=uigetfile('*.*');
img=imread(fullfile(img_path, img_name));


%% get the top-left and bottom-right corner of the rectangle window where
%% we perform TILT.
% by default, the first point should be top-left one and the second should
% be bottom-right. Don't mess up the order......
figure(1);
imshow(img);
hold on;
initial_points(:, 1) = ginput(1)';
plot(initial_points(1, 1), initial_points(2, 1), 'rx');
hold on; 
initial_points(:, 2) = ginput(1)';
plot(initial_points(1, 2), initial_points(2, 2), 'rx');
hold on; 

%% Run TILT with LADMAP and Warm Start Techniques
fprintf('\n');
[Dotau, A, E, f, tfm_matrix, focus_size, error_sign, UData, VData, XData, YData, A_scale]=...
    TILT(img, 'HOMOGRAPHY', initial_points, 'SAVE_PATH', fullfile(cd, task_name),...
    'BLUR_SIGMA', 2, 'BLUR_NEIGHBOR', 2, 'BRANCH', 0, 'NO_TRANSLATION', 0,...
    'PYRAMID_MAXLEVEL', 2, 'DISPLAY_INTER', 0, 'VWS', 1, 'SVDWS', 1);

%% Display the result
tfm=fliptform(maketform('projective', tfm_matrix'));
Dotau=imtransform(img, tfm, 'bilinear', 'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'Size', focus_size);
figure(3);
imshow(Dotau);
