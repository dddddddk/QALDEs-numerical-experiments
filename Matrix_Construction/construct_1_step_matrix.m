function L = construct_1_step_matrix(A, d)
n = length(A); k = length(d);
L = sparse(2*k*n, 2*k*n);
L(1:k*n, 1:k*n) = construct_forward_matrix(A, d, 'partial');
L(k*n+1:2*k*n, k*n+1:2*k*n) = construct_backward_matrix(A, d);
L(k*n+1:2*k*n, 1:k*n) = construct_T(n, k, '1');
end