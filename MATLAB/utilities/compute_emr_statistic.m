function T_O = compute_emr_statistic(y)
[m, n] = size(y);
yhat = [real(y)', imag(y)']';
Ryhat = 1/n * (yhat * yhat');

upperIdx = triu(true(2*m,2*m),1);
T_O = 1 + (1/m)*sum( Ryhat(upperIdx).^2 );

end