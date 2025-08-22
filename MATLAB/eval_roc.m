clearvars; clc;

addpath('constants');
addpath('functions');


p = 2;
m = 4;
n = 128;
SNR = -1;
H = get_H();
sigma1 = [1, 1, 1, 1];
sigma2 = [0.4, 0.8, 1.2, 1.6];
thresholds = 0:80;

[rao_Pd1, rao_Pfa1] = compute_one_bit_rao_roc(p, m, n, H, sigma1, SNR, thresholds, 1e6);
[rao_Pd2, rao_Pfa2] = compute_one_bit_rao_roc(p, m, n, H, sigma2, SNR, thresholds, 1e6);

[emr_Pd1, emr_Pfa1] = compute_one_bit_emr_roc(p, m, n, H, sigma1, SNR, thresholds, 1e6);
[emr_Pd2, emr_Pfa2] = compute_one_bit_emr_roc(p, m, n, H, sigma2, SNR, thresholds, 1e6);

figure
plot(rao_Pfa1, rao_Pd1, 'r', 'DisplayName', 'One-bit Rao, calibrated', 'LineWidth', 1.5); hold on
plot(rao_Pfa2, rao_Pd2, '--b', 'DisplayName', 'One-bit Rao, uncalibrated', 'LineWidth', 1.5);

plot(emr_Pfa1, emr_Pd1, 'k', 'DisplayName', 'One-bit EMR, calibrated', 'LineWidth', 1.5); hold on;
plot(emr_Pfa2, emr_Pd2, '--g', 'DisplayName', 'One-bit EMR, uncalibrated', 'LineWidth', 1.5);

grid on;
xlabel('P_{fa}');
ylabel('P_{d}');
legend('Location','northeast')
