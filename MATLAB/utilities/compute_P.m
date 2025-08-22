function P = compute_P(Rw, Rs, H)
    Rx = H * Rs * H' + Rw;
    D = diag(1 ./ sqrt(diag(Rx)));
    Px = D * Rx * D;
    P = [real(Px), -imag(Px); imag(Px), real(Px)];
end
