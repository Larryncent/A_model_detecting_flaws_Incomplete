function [L, S] = RobustPCA_alm(X, lambda, mu, tol, max_iter)
    % - method - algorithm method, given 'ALM', or 'only_penalty'
    % - X is a data matrix (of the size N x M) to be decomposed
    %   X can also contain NaN's for unobserved values
    % - lambda - regularization parameter, default = 1/sqrt(max(N,M))
    % - mu - the augmented lagrangian parameter, default = M*N/(4*normX)
    % - tol - reconstruction error tolerance, default = 1e-4
    % - max_iter - maximum number of iterations, default = 1000

    [M, N] = size(X);
    unobserved = isnan(X);
    X(unobserved) = 0;
    normX = norm(X, 'fro');

    % default arguments
    if nargin < 3
        lambda = 1 / sqrt(max(M,N));
    end
    if nargin < 4
        mu = M*N/(4*normX);
    end
    if nargin < 5
        tol = 1e-4;
    end
    if nargin < 6
        max_iter = 2000;
    end
    
    % initial solution
    L = zeros(M, N);
    S = zeros(M, N);
    Y = zeros(M, N);
    
    for iter = (1:max_iter)
        %ALM-----------------------------------------------------------
        L = SVT(1/mu, X - S + (1/mu)*Y);
        S = shrinkage_operator(lambda/mu, X - L + (1/mu)*Y);
        % and augmented lagrangian multiplier
        Z = X - L - S;
        Z(unobserved) = 0; % skip missing values
        Y = Y + mu*Z;
        %--------------------------------------------------------------
        err = norm(Z, 'fro') / normX;
        if (iter == 1) || (mod(iter, 10) == 0) || (err < tol)
            fprintf(1, 'iter: %04d\terr: %f\trank(L): %d\tcard(S): %d\n', ...
                    iter, err, rank(L), nnz(S(~unobserved)));
        end
        if (err < tol) break; end
    end
end

