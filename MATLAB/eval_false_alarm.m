clearvars; close all; clc;
addpath('functions');
addpath('utilities');

m = 4;
sigma = [1, 1, 1, 1];
thresholds = 0:80;
num_trials = 1e6;

false_alarm_analysis(m, [16, 32, 64], sigma, thresholds, num_trials);
