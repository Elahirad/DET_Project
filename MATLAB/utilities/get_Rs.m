function Rs = get_Rs(p, noise_power, SNR)
    Rs = eye(p) * (10 ^ (SNR / 10) * noise_power);
end
