function probabilities = compute_one_bit_emr_false_alarm_prob(m, n, sigmas, thresholds, iterations)
T_O = zeros(iterations, 1);
probabilities = zeros(size(thresholds));
for i=1:iterations
    x = zeros(m, n);

    for k=1:m
        x(k, :) = sqrt(sigmas(k)/2) * (randn(1,n) + 1i*randn(1,n));
    end

    y = sign(real(x)) + 1i*sign(imag(x));

    T_O(i) = m*n*(compute_emr_statistic(y)-1);
end
for i=1:size(thresholds, 2)
    probabilities(i) = mean(T_O > thresholds(i));
end
end