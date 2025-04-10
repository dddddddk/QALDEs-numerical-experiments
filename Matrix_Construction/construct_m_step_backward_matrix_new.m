function L = construct_m_step_backward_matrix_new(A, d, m)
n = length(A); k = length(d);
L = sparse((m*k + 1)*n, (m*k+1)*n);
W = construct_backward_matrix(A, d);
T2 = construct_T(n, k, '2'); 
W(1:n, :) = W(1:n, :) / sqrt(k);
L(1:n*k, 1:n*k) = W;
for i = 1:m-1
    L(n*i*k+1:n*(i+1)*k, n*i*k+1:n*(i+1)*k) = W;
    L(n*i*k+1:n*(i+1)*k, n*(i-1)*k+1:n*i*k) = T2 / sqrt(k);
end
L(n*m*k+1:n*m*k+n, n*(m-1)*k+1:n*m*k) = T2(1:n, :) / sqrt(k);
L(n*m*k+1:n*m*k+n, n*m*k+1:n*m*k+n) = eye(n) / sqrt(k);
end