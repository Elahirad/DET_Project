function false_alarm_analysis(m, n_arr, sigma_arr, thresholds, num_trials)

    figure

    for i = 1:length(n_arr)
        n = n_arr(i);
        rao_prob = compute_one_bit_rao_false_alarm_prob(m, n, sigma_arr, thresholds, num_trials);
        [alpha, beta] = get_beta_parameters(m, n);
        prob_beta = 1 - betacdf(get_normalized_thresholds(thresholds, m, n), alpha, beta);
        [k, ~] = get_chi_parameters(m, n);
        prob_chi = 1 - chi2cdf(thresholds, k);

        if i == 1
            semilogy(get_normalized_thresholds(thresholds, m, n), rao_prob, 'r', 'DisplayName', 'Simulation', 'LineWidth', 1.5); hold on
            semilogy(get_normalized_thresholds(thresholds, m, n), prob_beta, ':g', 'DisplayName', 'Beta approximation', 'LineWidth', 1.5);
            semilogy(get_normalized_thresholds(thresholds, m, n), prob_chi, '--b', 'DisplayName', 'Chi 2 approximation', 'LineWidth', 1.5);
        else
            semilogy(get_normalized_thresholds(thresholds, m, n), rao_prob, 'r', 'HandleVisibility', 'off', 'LineWidth', 1.5);
            semilogy(get_normalized_thresholds(thresholds, m, n), prob_beta, ':g', 'HandleVisibility', 'off', 'LineWidth', 1.5);
            semilogy(get_normalized_thresholds(thresholds, m, n), prob_chi, '--b', 'HandleVisibility', 'off', 'LineWidth', 1.5);
        end

    end

    grid on;
    xlabel('Threshold');
    ylabel('P_{fa}');
    ylim([1e-4, 1]);
    legend('Location', 'northeast')
end
