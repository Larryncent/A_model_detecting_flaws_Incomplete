% Modified by Xian-Hong Wang June 2023.
% Kelvin Zhang, Arvind Ganesh, February 2011. 
% Questions? zhangzdfaint@gmail.com, abalasu2@illinois.edu
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

function tau=tfm2para(tfm_matrix, XData, YData, mode)
% tfm2para will transpose tfm_matrix to its corresponding parameter.
% -------------------------input------------------------------------------
% tfm_matrix:       3-by-3 matrix.
% mode:             one of 'euclidean', 'euclidean_notranslation', 'affine', 'affine_notranslation', 
%                   'homography', 'homography_notranslation'
% -------------------------output-----------------------------------------
% tau:              p-by-1 real vector.
switch lower(mode)
    case 'affine'
        tau=zeros(6, 1);
        tau(1:3)=tfm_matrix(1, :)';
        tau(4:6)=tfm_matrix(2, :)';
    case {'homography'}
        tau=zeros(8, 1);
        tau(1)=tfm_matrix(2, 1)/tfm_matrix(1, 1);
        tau(2)=tfm_matrix(3, 1)/tfm_matrix(1, 1);
        tau(3)=tfm_matrix(1, 2)/tfm_matrix(2, 2);
        tau(4)=tfm_matrix(3, 2)/tfm_matrix(2, 2);
        pt=[XData; YData; ones(1, 2)];
        pt_trans=tfm_matrix*pt;
        pt_trans=pt_trans./(ones(3, 1)*pt_trans(3, :));
        pt_trans=pt_trans(1:2, :);
        tau(5:8)=pt_trans(:);
end