function T_L = compute_LMPIT_statistic(x)
    [m, n] = size(x);
    Rx = x*x'/n;
    T_L = real(n * trace((Rx*inv(diag(diag(Rx)))-eye(m))^2));
end