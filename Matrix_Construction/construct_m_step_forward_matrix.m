function L = construct_m_step_forward_matrix(A, d, m)
n = length(A); k = length(d);
% L = zeros(n*(k*m+1));
L = sparse(n*(k*m+1), n*(k*m+1));
M = construct_forward_matrix(A, d, 'partial');
T1 = construct_T(n, k, '1');
for i = 1:m-1
    L((i-1)*k*n+1:i*k*n, (i-1)*k*n+1:i*k*n) = M;
    L(i*k*n+1:(i+1)*k*n, (i-1)*k*n+1:i*k*n) = T1;
end
L((m-1)*k*n+1:(m*k+1)*n, (m-1)*k*n+1:(m*k+1)*n) = construct_forward_matrix(A, d, 'whole');
end