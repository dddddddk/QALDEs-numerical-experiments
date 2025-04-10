function A = construct_test_matrix_4(n)
A = sparse(n, n);
for k = 1:n-1
    A(k, k) = -0.01 * k^2;
    A(k, k+1) = 1;
end
A(n, n) = -0.01 * n^2;
end