function [Pd, Pfa] = compute_one_bit_rao_roc(p, m, n, H, sigmas, SNR, thresholds, iterations)
T_Rx = zeros(iterations, 1);
T_Rn = zeros(iterations, 1);
noise_power = mean(sigmas);
signal_sigma = sqrt(10^(SNR / 10) * noise_power);
Pd = zeros(size(thresholds));
Pfa = zeros(size(thresholds));
for i=1:iterations
    noise = diag(sqrt(sigmas)) * 1/sqrt(2) * (randn(m, n) + 1i*randn(m, n));

    s = (signal_sigma/sqrt(2)) * (randn(p,n) + 1i*randn(p,n));

    x = diag(sqrt(sigmas)) * H*s + noise;

    yx = sign(real(x)) + 1i*sign(imag(x));
    yn = sign(real(noise)) + 1i*sign(imag(noise));

    T_Rx(i) = compute_rao_statistic(yx);
    T_Rn(i) = compute_rao_statistic(yn);
end
for i=1:size(thresholds, 2)
    Pd(i) = mean(T_Rx > thresholds(i));
    Pfa(i) = mean(T_Rn > thresholds(i));
end
end