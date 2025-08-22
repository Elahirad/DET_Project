function thresholds = get_normalized_thresholds(thresholds, m, n)
    thresholds = thresholds / (m * n * (m - 1));
end
