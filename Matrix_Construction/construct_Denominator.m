function D = construct_Denominator(A, d)
% A: a square matrix
% d: a coefficients list
n = length(A);
I = eye(n);
k = length(d);
D = d(k) * I;
for i = 1:k-1
    D = A * D + d(k-i) * I;
end
end