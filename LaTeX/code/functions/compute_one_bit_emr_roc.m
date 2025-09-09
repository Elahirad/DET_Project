function [Pd, Pfa] = compute_one_bit_emr_roc(p, m, n, H, sigmas, SNR, thresholds, iterations)
T_Ox = zeros(iterations, 1);
T_On = zeros(iterations, 1);
noise_power = sum(sigmas) / m;
signal_sigma = sqrt(10^(SNR / 10) * noise_power);
Pd = zeros(size(thresholds));
Pfa = zeros(size(thresholds));
for i=1:iterations
    noise = diag(sqrt(sigmas)) * 1/sqrt(2) * (randn(m, n) + 1i*randn(m, n));

    s = (signal_sigma/sqrt(2)) * (randn(p,n) + 1i*randn(p,n));

    x = diag(sqrt(sigmas)) * H*s + noise;

    yx = sign(real(x)) + 1i*sign(imag(x));
    yn = sign(real(noise)) + 1i*sign(imag(noise));

    T_Ox(i) = m*n*(compute_emr_statistic(yx) - 1);
    T_On(i) = m*n*(compute_emr_statistic(yn) - 1);
end
for i=1:size(thresholds, 2)
    Pd(i) = mean(T_Ox > thresholds(i));
    Pfa(i) = mean(T_On > thresholds(i));
end
end