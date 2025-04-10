function d = Taylor_Coefficient(k)
% This function return the coefficients of taylor approximants to exponential
% function. The polynomial order is k.
d = zeros(k+1, 1);
d(1) = 1;
for i = 2:k+1
    d(i) = d(i-1) / (i-1);
end
end