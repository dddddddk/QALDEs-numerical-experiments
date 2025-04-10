function L = construct_m_step_forbackward_matrix(A, d, m)
n = length(A); k = length(d);
L = sparse((m+1)*k*n, (m+1)*k*n);
L(1:n*k, 1:n*k) = construct_forward_matrix(A, d, 'partial');
L(n*k+1:2*n*k, 1:n*k) = construct_T(n, k, '1');
W = construct_backward_matrix(A, d);
T2 = construct_T(n, k, '2');
L(n*k+1:2*n*k, n*k+1:2*n*k) = W;
for i = 2:m
    L(n*i*k+1:n*(i+1)*k, n*i*k+1:n*(i+1)*k) = W;
    L(n*i*k+1:n*(i+1)*k, n*(i-1)*k+1:n*i*k) = T2;
end
end