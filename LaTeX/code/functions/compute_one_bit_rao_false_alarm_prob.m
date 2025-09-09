function probabilities = compute_one_bit_rao_false_alarm_prob(m, n, sigmas, thresholds, iterations)
T_R = zeros(iterations, 1);
probabilities = zeros(size(thresholds));
for i=1:iterations
    x = zeros(m, n);

    for k=1:m
        x(k, :) = sqrt(sigmas(k)/2) * (randn(1,n) + 1i*randn(1,n));
    end

    y = sign(real(x)) + 1i*sign(imag(x));

    T_R(i) = compute_rao_statistic(y);
end
for i=1:size(thresholds, 2)
    probabilities(i) = mean(T_R > thresholds(i));
end
end