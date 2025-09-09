function [Pd_L, Pd_O, Pd_R, Pd_R2_47, Pd_E] = compute_detection_prob_vs_SNR(p, m, n, sigmas, SNR_arr, target_Pfa, iterations)

Pd_L=zeros(size(SNR_arr));
Pd_O=zeros(size(SNR_arr));
Pd_R=zeros(size(SNR_arr));
Pd_R2_47=zeros(size(SNR_arr));
Pd_E=zeros(size(SNR_arr));

T_En = zeros(iterations, 1);
parfor i=1:iterations
    noise = diag(sqrt(sigmas)) * 1/sqrt(2) * (randn(m, n) + 1i*randn(m, n));
    T_En(i) = compute_inf_emr_statistic(noise);
end
thr_O = chi2inv(1-target_Pfa, 2*m*(m-1))/(m*n)+1;
thr_E = quantile(T_En, 1-target_Pfa);
[k, ~] = get_chi_parameters(m, n);

thr_lmp = chi2inv(1-target_Pfa, k);
thr_rao = thr_lmp;
thr_rao247 = thr_lmp;

noise_power = mean(sigmas);
D = diag(sqrt(sigmas));
snr_lin = 10.^(SNR_arr/10);
h = waitbar(0, 'Processing...');
for j=1:numel(SNR_arr)
    signal_sigma = sqrt(snr_lin(j) * noise_power);
    T_Ex = zeros(iterations, 1);
    T_Ox = zeros(iterations, 1);
    parfor i=1:iterations
        noise = D * 1/sqrt(2) * (randn(m, n) + 1i*randn(m, n));
        s = (signal_sigma/sqrt(2)) * (randn(p,n) + 1i*randn(p,n));


        H = (randn(m,p) + 1i*randn(m,p))/sqrt(2);
        H = H ./ vecnorm(H,2,1);

        x = D * H*s + noise;
        x_q = sign(real(x)) + 1i*sign(imag(x));


        T_Ex(i) = compute_inf_emr_statistic(x);
        T_Ox(i) = compute_emr_statistic(x_q);
    end
    Pd_O(j) = mean(T_Ox > thr_O);
    Pd_E(j) = mean(T_Ex > thr_E);
    waitbar(j / length(SNR_arr), h, sprintf('Processing... %d%%', round(j/length(SNR_arr)*100)));
end
close(h);

H = get_H();
for j=1:numel(SNR_arr)
    signal_sigma = sqrt(10^(SNR_arr(j) / 10) * noise_power);
    P = compute_P(diag(sigmas), eye(p)*signal_sigma^2, H);
    [k, delta] = get_chi_parameters(m, n, P);
    delta_lmp = pi^2 / 4 * delta;
    [k247, delta247] = get_chi_parameters(m, round(2.47*n), P);


    Pd_L(j) = 1 - ncx2cdf(thr_lmp, k, delta_lmp);
    Pd_R(j) = 1 - ncx2cdf(thr_rao, k, delta);
    Pd_R2_47(j) = 1 - ncx2cdf(thr_rao247, k247, delta247);
    
end
end