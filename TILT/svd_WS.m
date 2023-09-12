% Xiang Ren, Sep. 2011.
%
% Questions? renxiangzju@gmail.com
%
% Copyright: Microsoft Research Asia, Beijing
%
% Reference: Linearized Alternating Direction Method with Adaptive Penalty
%            for Fast Solving Transform Invariant Low-rank Texture.
%            Xiang Ren, Zhouchen Lin. Submitted to IJCV.
%

function [U, S, V] = svd_WS(W, U0, S0, V0)

% This matlab code implements svd using warm-start technique
%----------------------------------------------------------------------
% min 1/2||W - USV'||_F^2
% s.t., U'U=I, V'V=I
%--------------------------------
% inputs:
%        W -- M*N data matrix
%        U0, S0, V0 -- previous decomposed matrices
%        
% outputs:
%        U, S, V -- new decomposed matrices
%

%% Initilization (4 times multiplications)
display = 0;
[m n] = size(W);
U = U0;
V = V0;
S = S0;
USVT = U*S*V';
USVTWT = USVT*W';
A = USVTWT - USVTWT'; % m-by-m
VSUTW = USVT'*W;
B = VSUTW - VSUTW';% n-by-n
r = size(U,2);
% calculate diag(U'WV):
% diagUTWV = diag(diag(U'*W*V));
diagUTWV = zeros(r,r);
for i = 1: r
    diagUTWV(i,i) = U(:,i)'*W*V(:,i);
end
diagS = S-diagUTWV;
IA = eye(m);
IB = eye(n);

%% Compute the derivatives at zero (9 times multiplication)
h0 = W - USVT;
UDSVT = U*diagS*V';
g0 = A*USVT + UDSVT - USVT*B;
gg0 = -A*(g0+UDSVT) + (g0+UDSVT)*B;
f0 = norm(h0,'fro')^2/2;
f0_first = trace(h0'*g0);
f0_second = trace(g0'*g0 + h0'*gg0);

%% Compute optimal step-size
test = f0_first^2-2*f0_second*f0;
if f0_first^2-2*f0_second*f0 < 0
    t = - f0_first/f0_second; % - b/2a
else
    poly = [f0_second/2 f0_first f0];
    t_root = roots(poly);
    t = t_root(1);
end
% disp(['value of t: ' num2str(t)]);


%% Update the matrices  (4 times multiplications, 2 inverses)
invA = inv(IA+t*A/2);
invB = inv(IB+t*B/2);
U = invA*(IA-t*A/2)*U;
V = invB*(IB-t*B/2)*V;
S = S - t*diagS;

%% Compute the error
if display == 1
ht = W -U*S*V';
ft = norm(ht,'fro')^2/2;
disp(['ft: ' num2str(ft)]);
end

end

