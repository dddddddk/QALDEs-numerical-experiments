function L = construct_m_step_forbackward_matrix_new(A, d, m)
% This function construct a linear system different from the m step
% forbackward matrix. It is preconditioned by multiplying 1/sqrt(k) in the
% first row of T1, T2, and W.
n = length(A); k = length(d);
L = sparse((m+1)*k*n, (m+1)*k*n);
L(1:n*k, 1:n*k) = construct_forward_matrix(A, d, 'partial');
L(n*k+1:2*n*k, 1:n*k) = construct_T(n, k, '1') / sqrt(k);
W = construct_backward_matrix(A, d);
W(1:n, :) = W(1:n, :) / sqrt(k);
T2 = construct_T(n, k, '2') / sqrt(k);
L(n*k+1:2*n*k, n*k+1:2*n*k) = W;
for i = 2:m
    L(n*i*k+1:n*(i+1)*k, n*i*k+1:n*(i+1)*k) = W;
    L(n*i*k+1:n*(i+1)*k, n*(i-1)*k+1:n*i*k) = T2;
end