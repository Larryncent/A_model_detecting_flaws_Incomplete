close all;
clear all;
clc;
%% load the image
[img_name, img_path]=uigetfile('*.*');
img=imread(fullfile(img_path, img_name));
%img = imread('./test_images/building.png');
%img = imresize(img, [256, 256]);
[M, N, P] = size(img);
img = double(img);
img = img/255;
L = zeros(M, N,3);

%Input mu 
list_mu = [40, 20, 30];
[tmp1, numbers_mu] = size(list_mu); 

%Input lambda
list_lambda = [0.004, 0.002];
[tmp2, numbers_lambda] = size(list_lambda); 


for j = 1:numbers_mu
    for k = 1:numbers_lambda
        for i = 1:3
            fprintf('C %i \n', i);
            [l, S] = RobustPCA_alm(img(:, :, i), list_lambda(k), list_mu(j));
            L(:,:, i) = l;
            save_path = './results_alm/';
            filename = strcat('lda_', string(list_lambda(k)), ...
                'mu_', string(list_mu(j)), '.jpg');
            imwrite(L,strcat(save_path, filename))
        end
    end
end    




