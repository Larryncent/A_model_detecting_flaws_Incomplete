function r = SVT(tau, X)
    % Singular value 
    [U, S, V] = svd(X, 'econ');
    r = U*shrinkage_operator(tau, S)*V';
end