function r = shrinkage_operator(tau, X)
    % shrinkage_operator
    r = sign(X) .* max(abs(X) - tau, 0);
end