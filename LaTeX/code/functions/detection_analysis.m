function detection_analysis(Rs, Rw, n_arr, H, sigma_arr, SNR, thresholds, num_trials, x_range)
    p = size(Rs, 1);
    m = size(Rw, 1);

    P = compute_P(Rw, Rs, H);

    figure
    disp(["SNR=", num2str(SNR)]);
    for i = 1:length(n_arr)
        n = n_arr(i);
        rao_prob = compute_one_bit_rao_detection_prob(p, m, n, H, sigma_arr, SNR, thresholds, num_trials);
        [alpha, beta] = get_beta_parameters(m, n, P);
        prob_beta = 1 - betacdf(get_normalized_thresholds(thresholds, m, n), alpha, beta);
        [k, delta] = get_chi_parameters(m, n, P);
        prob_chi = 1 - ncx2cdf(thresholds, k, delta);

        if i == 1
            plot(get_normalized_thresholds(thresholds, m, n), rao_prob, 'r', 'DisplayName', 'Simulation', 'LineWidth', 1.5); hold on
            plot(get_normalized_thresholds(thresholds, m, n), prob_beta, ':g', 'DisplayName', 'Beta approximation', 'LineWidth', 1.5);
            plot(get_normalized_thresholds(thresholds, m, n), prob_chi, '--b', 'DisplayName', 'Chi 2 approximation', 'LineWidth', 1.5);
        else
            plot(get_normalized_thresholds(thresholds, m, n), rao_prob, 'r', 'HandleVisibility', 'off', 'LineWidth', 1.5);
            plot(get_normalized_thresholds(thresholds, m, n), prob_beta, ':g', 'HandleVisibility', 'off', 'LineWidth', 1.5);
            plot(get_normalized_thresholds(thresholds, m, n), prob_chi, '--b', 'HandleVisibility', 'off', 'LineWidth', 1.5);
        end

        disp(['Beta distribution error @ n=', num2str(n), '->', num2str(mean((prob_beta-rao_prob).^2))]);
        disp(['Chi-2 distribution error@ n=', num2str(n), '->', num2str(mean((prob_chi-rao_prob).^2))]);

    end

    grid on;
    xlim(x_range);
    xlabel('Threshold');
    ylabel('P_{d}');
    legend('Location', 'northeast')
end
