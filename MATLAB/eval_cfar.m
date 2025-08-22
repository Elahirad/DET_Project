clearvars; close all; clc;
addpath('functions');

m = 4;
n = 128;
sigma1 = [1, 1, 1, 1];
sigma2 = [0.4, 0.8, 1.2, 1.6];
sigma3 = [0.5, 0.75, 1.2, 1.5];
thresholds = 0:80;

rao_prob1 = compute_one_bit_rao_false_alarm_prob(m, n, sigma1, thresholds, 1e6);
rao_prob2 = compute_one_bit_rao_false_alarm_prob(m, n, sigma2, thresholds, 1e6);
rao_prob3 = compute_one_bit_rao_false_alarm_prob(m, n, sigma3, thresholds, 1e6);

emr_prob1 = compute_one_bit_emr_false_alarm_prob(m, n, sigma1, thresholds, 1e6);
emr_prob2 = compute_one_bit_emr_false_alarm_prob(m, n, sigma2, thresholds, 1e6);
emr_prob3 = compute_one_bit_emr_false_alarm_prob(m, n, sigma3, thresholds, 1e6);

figure
semilogy(thresholds, rao_prob1, '-.sk', 'DisplayName', 'One-bit Rao, case 1'); hold on
semilogy(thresholds, rao_prob2, '-.m^', 'DisplayName', 'One-bit Rao, case 2');
semilogy(thresholds, rao_prob3, '-.yo', 'DisplayName', 'One-bit Rao, case 3');

semilogy(thresholds, emr_prob1, '-sr', 'DisplayName', 'One-bit EMR, case 1');
semilogy(thresholds, emr_prob2, '-^g', 'DisplayName', 'One-bit EMR, case 2');
semilogy(thresholds, emr_prob3, '-ob', 'DisplayName', 'One-bit EMR, case 3');

ylim([1e-4, 1]);
grid on;
xlabel('Threshold');
ylabel('P_{fa}');
legend('Location','northeast')
