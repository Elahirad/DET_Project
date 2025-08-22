function T_EMR = compute_inf_emr_statistic(x)
[m, n] = size(x);
Rx = (x*x')/n;
num = norm(Rx,'fro')^2;
den = (trace(Rx)/m)^2;
T_EMR = (num/m) / den;
end