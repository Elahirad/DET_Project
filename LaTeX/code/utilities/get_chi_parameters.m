function [k, delta] = get_chi_parameters(m, n, P)

    if nargin < 3
        k = m ^ 2 - m;
        delta = 0;
    else
        Px = P(1:m, 1:m) + 1i*P(m+1:2*m, 1:m);
        k = m ^ 2 - m;
        delta = 4 / pi ^ 2 * n * real(trace((Px - eye(m)) * (Px - eye(m)')));
    end

end
