clearvars; close all; clc;
addpath('functions');
addpath('utilities');

%% SNR = -9dB
p = 2;
m = 4;
SNR = -9;
H = get_H();
sigma = [1, 1, 1, 1];
thresholds = 0:80;
num_trials = 1e6;
Rw = get_Rw(sigma);
Rs = get_Rs(p, sum(sigma .^ 2) / m, SNR);

detection_analysis(Rs, Rw, [64, 128, 256], H, sigma, SNR, thresholds, num_trials, [0, 0.04]);

%% SNR = 4dB

p = 2;
m = 4;
SNR = 4;
H = get_H();

sigma = [1, 1, 1, 1];
thresholds = 0:240;
num_trials = 1e6;
Rw = get_Rw(sigma);
Rs = get_Rs(p, sum(sigma .^ 2) / m, SNR);

detection_analysis(Rs, Rw, [64, 128, 256], H, sigma, SNR, thresholds, num_trials, [0, 0.12]);
