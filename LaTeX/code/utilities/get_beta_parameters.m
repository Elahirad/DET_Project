function [alpha, beta] = get_beta_parameters(m, n, P)

    if nargin < 3
        alpha = (n * m * (m - 1) - 2) / (2 * n);
        beta = ((n - 1) * (n * m * (m - 1) - 2)) / (2 * n);
    else

        iprime = @(i) (i <= m) .* (i + m) + (i > m) .* (i - m);
        clip = @(A) min(max(A, -1), 1);

        P = clip(P);
        h = (2 / pi) * asin(P);

        An2 = n * (n - 1);
        mu1 = 0;

        for i = 1:m

            for j = i + 1:m
                gij = (4 / n ^ 2) * (n + An2 * (h(i, j) ^ 2 + h(iprime(i), j) ^ 2));
                mu1 = mu1 + gij;
            end

        end

        mu1 = mu1 / (2 * m * (m - 1));

        habcd = @(a, b, c, d) ...
            (16 * mvncdf(zeros(1, 4), inf(1, 4), zeros(1, 4), P([a b c d], [a b c d])) ...
            - 1 ...
            - (h(a, b) + h(a, c) + h(a, d) + h(b, c) + h(b, d) + h(c, d)));

        v1 = @(i, j, k, l) (habcd(i, j, k, l) + habcd(i, j, iprime(k), iprime(l)) ...
            + habcd(iprime(i), iprime(j), k, l) + habcd(iprime(i), iprime(j), iprime(k), iprime(l)));
        v2 = @(i, j, k, l) (habcd(iprime(i), j, iprime(k), l) - habcd(iprime(i), j, k, iprime(l)) ...
            - habcd(i, iprime(j), iprime(k), l) + habcd(i, iprime(j), k, iprime(l)));
        v3 = @(i, j, k, l) (habcd(iprime(i), j, k, l) + habcd(iprime(i), j, iprime(k), iprime(l)) ...
            - habcd(i, iprime(j), k, l) - habcd(i, iprime(j), iprime(k), iprime(l)));
        v4 = @(i, j, k, l) (habcd(i, j, iprime(k), l) - habcd(i, j, k, iprime(l)) ...
            + habcd(iprime(i), iprime(j), iprime(k), l) - habcd(iprime(i), iprime(j), k, iprime(l)));

        An2 = n * (n - 1);
        An3 = n * (n - 1) * (n - 2);

        tau1 = @(i, j) ...
            16 / n ^ 4 * An2 * (1 + habcd(i, iprime(i), j, iprime(j)) ^ 2) ...
            + 32 / n ^ 4 * An3 * ((h(i, j) ^ 2 + h(iprime(i), j) ^ 2) ...
            + habcd(i, iprime(i), j, iprime(j)) * (h(i, j) ^ 2 - h(iprime(i), j) ^ 2));

        tau2 = @(i, j, k) ...
            (4 / n ^ 4) * An2 * (4 * (h(j, k) ^ 2 + h(j, iprime(k)) ^ 2) ...
            + (habcd(i, iprime(i), j, iprime(k)) + habcd(i, iprime(i), iprime(j), k)) .^ 2 ...
            + (habcd(i, iprime(i), j, k) - habcd(i, iprime(i), iprime(j), iprime(k))) .^ 2) ...
            + (16 / n ^ 4) * An3 * ((habcd(i, iprime(i), j, k) - habcd(i, iprime(i), iprime(j), iprime(k))) ...
            .* (h(i, k) .* h(iprime(i), j) + h(i, j) .* h(iprime(i), k)) ...
            + (habcd(i, iprime(i), j, iprime(k)) + habcd(i, iprime(i), iprime(j), k)) ...
            .* (h(i, k) .* h(i, j) - h(iprime(i), j) .* h(iprime(i), k)) ...
            + 2 * h(j, k) .* (h(i, k) .* h(i, j) + h(iprime(i), j) .* h(iprime(i), k)) ...
            + 2 * h(j, iprime(k)) .* (h(i, k) .* h(iprime(i), j) - h(i, j) .* h(iprime(i), k)));

        tau3 = @(i, j, k, l) ...
            (2 / n ^ 4) * An2 * (v1(i, j, k, l) .^ 2 + v2(i, j, k, l) .^ 2 + v3(i, j, k, l) .^ 2 + v4(i, j, k, l) .^ 2) ...
            + (16 / n ^ 4) * An3 * (v1(i, j, k, l) .* h(i, j) .* h(k, l) ...
            + v2(i, j, k, l) .* h(iprime(i), j) .* h(iprime(k), l) ...
            + v3(i, j, k, l) .* h(k, l) .* h(iprime(i), j) ...
            + v4(i, j, k, l) .* h(i, j) .* h(k, iprime(l)));

        sig1 = 0;

        for i = 1:m

            for j = i + 1:m

                for k = 1:m

                    for l = k + 1:m

                        gijkl = (32 * (n - 1) * (2 * n - 3) / n ^ 3) * ...
                            (h(i, j) ^ 2 + h(iprime(i), j) ^ 2) * (h(k, l) ^ 2 + h(iprime(k), l) ^ 2);

                        if (i == k && j == l)
                            fijkl = tau1(i, j);
                        elseif (i == k && j ~= l)
                            fijkl = tau2(i, j, l);
                        elseif (i == l)
                            fijkl = tau2(i, j, k);
                        elseif (j == k)
                            fijkl = tau2(j, i, l);
                        elseif (j == l && i ~= k)
                            fijkl = tau2(j, i, k);
                        else
                            fijkl = tau3(i, j, k, l);
                        end

                        sig1 = sig1 + (fijkl - gijkl);
                    end

                end

            end

        end

        sigma1_sq = sig1 / (4 * m ^ 2 * (m - 1) ^ 2);

        alpha = (mu1 * (mu1 - mu1 ^ 2 - sigma1_sq)) / sigma1_sq;
        beta = ((1 - mu1) * (mu1 - mu1 ^ 2 - sigma1_sq)) / sigma1_sq;
    end
