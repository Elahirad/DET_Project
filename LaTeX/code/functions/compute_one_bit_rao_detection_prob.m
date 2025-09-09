function probabilities = compute_one_bit_rao_detection_prob(p, m, n, H, sigmas, SNR, thresholds, iterations)
T_R = zeros(iterations, 1);
noise_power = sum(sigmas.^2) / m;
signal_sigma = sqrt(10^(SNR / 10) * noise_power);
probabilities = zeros(size(thresholds));
for i=1:iterations
    noise = zeros(m, n);


    for k=1:m
        noise(k, :) = sqrt(sigmas(k)/2) * (randn(1,n) + 1i*randn(1,n));
    end
    s = (signal_sigma/sqrt(2)) * (randn(p,n) + 1i*randn(p,n));

    x = H*s + noise;

    y = sign(real(x)) + 1i*sign(imag(x));

    T_R(i) = compute_rao_statistic(y);
end
for i=1:size(thresholds, 2)
    probabilities(i) = mean(T_R > thresholds(i));
end
end