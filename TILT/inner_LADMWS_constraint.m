% Xiang Ren, Aug. 2011.
%
% Questions? renxiangzju@gmail.com
%
% Copyright: Microsoft Research Asia, Beijing
%
% Reference: Linearized Alternating Direction Method with Adaptive Penalty
%            for Fast Solving Transform Invariant Low-rank Texture.
%            Xiang Ren, Zhouchen Lin. Submitted to IJCV.
%
function [A, E, delta_tau, Y, mu, f, error_sign] = ...
    inner_LADMWS_constraint(Dotau, A0, E0, Y0, mu0, J, J_S, inner_para)

% This matlab code implements Linearized ADM with Adaptive Penalty (LADMAP+WS) method for TILT inner loop
%----------------------------------------------------------------------
% min |A|_*+lambda*|E|_1
% s.t., A + E = Do\tau + J\Delta\tau
%--------------------------------
% inputs:
%        Do\tau -- M*N data matrix
%        J -- M*N*P tensor, Jacobi
%        J_S -- linear constriant on \Delta\tau, k-by-p matrix,
%               k is number of constriants.
%        tol1 -- first threshold parameter for inner loop constraint
%        tol2 -- second threshold parameter for inner loop constraint
%        svd_tol -- threshold for SVD warm start
%
% outputs:
%        A -- M*N low-rank matrix
%        E -- M*N sparse error matrix
%
% ---------------------------------------------------------------------

%% parse parameters
[m n] = size(Dotau);
if (isempty(inner_para.inner_disp))
    inner_disp = 100; 
else
    inner_disp = inner_para.inner_disp;
end
if (isempty(inner_para.tol1))
    tol1 = 1e-7; 
else
    tol1 = inner_para.tol1;
end
if (isempty(inner_para.tol2))
    tol2 = 1e-2; 
else
    tol2 = inner_para.tol2;
end
if (isempty(inner_para.mu_max))
    mu_max = 1e10; 
else
    mu_max = inner_para.mu_max;
end
if (isempty(inner_para.rho))
    rho = 8; 
else
    rho = inner_para.rho;
end
c = 1;
if (isempty(inner_para.lambda))
    lambda = c/sqrt(m);
else
    lambda = inner_para.lambda;
end
if (isempty(inner_para.gamma))
    gamma = 1; 
else
    gamma = inner_para.gamma;
end

if (isempty(inner_para.svdws))
    svdws = 0; 
else
    svdws = inner_para.svdws;
end

% Whether do svd_WS
if svdws == 1
    svd_tol = 3e-3; % can be tuned
else
    svd_tol = -1;
end

%% Initilization
invJS = inv(J'*J+J_S'*J_S); % p-by-p matrix
Dotau = Dotau(:); % turn to vector
A = A0;
E = E0;
Y = Y0; % Lagrange multiplier
mu = mu0;
W = 0; % can be A0 I think!!
Normalization = norm(Dotau, 2);

%% Begin inner loop
inner = 0;
constriant = J*(invJS*(J'*(A+E-Dotau)));
constriant = A+E-Dotau-constriant; %initial value of "J^\perp(A+E-Dotau)"
while 1
    inner = inner+1;
    % update E:
    E_old = E;
    temp = E-Y/mu-constriant;
    E = sign(temp).*((abs(temp)-lambda/mu).*(abs(temp)>lambda/mu));
    % update A:
    A_old = A;
    temp = J*(invJS*(J'*(A+E-Dotau))); % computes (2kpmn+kp^2) times
    temp = A+E-Dotau-temp;
    temp = A-Y/mu-temp;
    W_old = W;
    W = reshape(temp,m,n); % reshape b ack to m-by-n to do SVD
    % do SVT
    if norm(W-W_old,'fro')/norm(W_old,'fro')> svd_tol
        [U S V] = svd(W,'econ');
    else
        [U S V] = svd_WS(W,U,S,V);
    end
    A = U*((S>1/mu).*(S-1/mu))*V';
    A = A(:); 
    % compute diff
    constriant = J*(invJS*(J'*(A+E-Dotau))); 
    constriant = A+E-Dotau-constriant;
    relE = norm(E-E_old, 2)/Normalization;
    relA = norm(A-A_old, 2)/Normalization;
    Err = norm(constriant,2)/Normalization; % new constriant error
    % display results
    if mod(inner,inner_disp)==0 
        disp(['iter ' num2str(inner) ',mu=' num2str(mu) ...
             ',||A||_*=' num2str(sum(sum(abs((S>1/mu).*(S-1/mu))))) ...
             ',consErr=' num2str(Err)]);
    end
    % update Y:
    Y = Y + gamma*mu*(constriant); 
    % stop criterion
    if Err <tol1 % && max(relA,relE) < tol2
        break;
    end
    % update mu:
    if min(mu,1)*max(relA,relE) < tol2
            mu = min(mu_max, mu*rho);
            TEST_MU(inner) = mu;     
    end
end
% compute \Delta\tau
M = A+E-Dotau;
delta_tau = invJS*(J'*M);
% computer objective function value
f = sum(sum(abs((S>1/mu).*(S-1/mu))))+lambda*sum(sum(abs(E)));
error_sign = 0;



