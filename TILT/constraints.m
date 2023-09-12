% Modified by Xian-Hong Wang June 2023.
% Kelvin Zhang, Arvind Ganesh, February 2011. 
%
% Reference: TILT: Transform Invariant Low-rank Textures  
%            Zhengdong Zhang, Xiao Liang, Arvind Ganesh, and Yi Ma. Proc. of ACCV, 2010.
%

function S=constraints(tau, XData, YData, mode)
% constraints() will get the linearize constraints of tau according to
% mode.
% -----------------------------input--------------------------------------
% tau:          p-by-1 real vector.
% mode:         one of 'affine', 'affine_notranslation', 'homography',
%              
% ----------------------------output--------------------------------------
% S:            linearized constraints on tau.
switch lower(mode)
    case 'affine'
        S=zeros(2, 6);
        vec1=[tau(1); tau(4)];
        vec2=[tau(2); tau(5)];
        V=sqrt(norm(vec1)^2*norm(vec2)^2-(vec1'*vec2)^2);
        S(1, 1)=(tau(1)*norm(vec2)^2-vec1'*vec2*tau(2))/V;
        S(1, 2)=(tau(2)*norm(vec1)^2-vec1'*vec2*tau(1))/V;
        S(1, 4)=(tau(4)*norm(vec2)^2-vec1'*vec2*tau(5))/V;
        S(1, 5)=(tau(5)*norm(vec1)^2-vec1'*vec2*tau(4))/V;
        S(2, 1)=2*tau(1);
        S(2, 2)=-2*tau(2);
        S(2, 4)=2*tau(4);
        S(2, 5)=-2*tau(5);
    case 'homography'
        S=zeros(0, 8);
end