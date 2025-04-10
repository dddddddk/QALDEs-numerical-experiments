function W = construct_backward_matrix(A, d)
% A: a square matrix
% d: a coefficients list
k = length(d); n = length(A);
W = sparse(k*n, k*n); I = eye(n);
W(1:n, 1:n) = I;
for i = 1:k-1
    W(1:n, n*i+1:n*(i+1)) = I;
    W(n*i+1:n*(i+1), n*i+1:n*(i+1)) = d(k+1-i)/d(k-i) * A;
    W(n*i+1:n*(i+1), n*(i-1)+1:n*i) = I;
end
end