function T_R = compute_rao_statistic(y)
[m, n] = size(y);
Ry = 1/n * (y * y');


upperIdx = triu(true(m,m),1);
T_R = n / 2 *sum( abs(Ry(upperIdx)).^2) ;

end