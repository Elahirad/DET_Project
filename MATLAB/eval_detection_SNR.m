%%
clearvars; clc;

addpath('functions');
addpath('utilities');
%%
p = 2;
m = 4;
n = 128;
sigma1 = [1, 1, 1, 1];

SNR_res = 15;
SNR_start = -7;
SNR_end = 7;
SNR_range = SNR_end - SNR_start + 1;

SNR = linspace(SNR_start, SNR_end, SNR_res);
[Pd_L, Pd_O, Pd_R, Pd_R2_47, Pd_E] = compute_detection_prob_vs_SNR(p, m, n, sigma1, SNR, 1e-4, 1e6);
Pd_E_shift = circshift(Pd_E, round(2* SNR_res/ SNR_range));
Pd_E_shift(1:round(2* SNR_res/ SNR_range)) = 0;
Pd_L_shift = circshift(Pd_L, round(2* SNR_res/ SNR_range));
Pd_L_shift(1:round(2* SNR_res/ SNR_range)) = 0;

figure
plot(SNR, Pd_L, 'r', 'DisplayName', 'LMPIT', 'LineWidth', 1.5); hold on;
plot(SNR, Pd_L_shift, '--b', 'DisplayName', 'LMPIT,2dB shift', 'LineWidth', 1.5);
plot(SNR, Pd_O, 'k', 'DisplayName', 'One-Bit EMR', 'LineWidth', 1.5);
plot(SNR, Pd_R, 'y', 'DisplayName', 'One-Bit Rao', 'LineWidth', 1.5);
plot(SNR, Pd_R2_47, '--c', 'DisplayName', 'One-Bit Rao,2.47n', 'LineWidth', 1.5);
plot(SNR, Pd_E, 'g', 'DisplayName', 'inf-EMR', 'LineWidth', 1.5);
plot(SNR, Pd_E_shift, '--m', 'DisplayName', 'inf-EMR,2dB shift', 'LineWidth', 1.5);


grid on;
xlabel('SNR');
ylabel('P_{d}');
legend('Location','northeast')

%%
p = 2;
m = 4;
n = 2048;
sigma1 = [1, 1, 1, 1];

SNR_res = 15;
SNR_start = -15;
SNR_end = -1;
SNR_range = SNR_end - SNR_start + 1;
SNR = linspace(SNR_start, SNR_end, SNR_res);


[Pd_L, Pd_O, Pd_R, Pd_R2_47, Pd_E] = compute_detection_prob_vs_SNR(p, m, n, sigma1, SNR, 1e-4, 1e6);
Pd_E_shift = circshift(Pd_E, round(2* SNR_res/ SNR_range));
Pd_E_shift(1:round(2* SNR_res/ SNR_range)) = 0;
Pd_L_shift = circshift(Pd_L, round(2* SNR_res/ SNR_range));
Pd_L_shift(1:round(2* SNR_res/ SNR_range)) = 0;

figure
plot(SNR, Pd_L, 'r', 'DisplayName', 'LMPIT', 'LineWidth', 1.5); hold on;
plot(SNR, Pd_L_shift, '--b', 'DisplayName', 'LMPIT,2dB shift', 'LineWidth', 1.5);
plot(SNR, Pd_O, 'k', 'DisplayName', 'One-Bit EMR', 'LineWidth', 1.5);
plot(SNR, Pd_R, 'y', 'DisplayName', 'One-Bit Rao', 'LineWidth', 1.5);
plot(SNR, Pd_R2_47, '--c', 'DisplayName', 'One-Bit Rao,2.47n', 'LineWidth', 1.5);
plot(SNR, Pd_E, 'g', 'DisplayName', 'inf-EMR', 'LineWidth', 1.5);
plot(SNR, Pd_E_shift, '--m', 'DisplayName', 'inf-EMR,2dB shift', 'LineWidth', 1.5);


grid on;
xlabel('SNR');
ylabel('P_{d}');
legend('Location','northeast')

