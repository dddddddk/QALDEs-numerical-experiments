function M = construct_forward_matrix(A, d, type)
% A: a square matrix
% d: a coefficients list
% type: if type == 'whole' then the function return the whole forward
% matrix, else the function return the forward matrix without the last row
% and last column.
k = length(d); n = length(A);
M = sparse((k+1)*n, (k+1)*n); I = eye(n);
M(1:n, 1:n) = I;
M(k*n+1:(k+1)*n, 1:n) = -I;
for i = 1:k-1
    M(i*n+1:(i+1)*n, i*n+1:(i+1)*n) = I;
    M(i*n+1:(i+1)*n, (i-1)*n+1:(i)*n) = -d(i+1)/d(i) * A;
    M(k*n+1:(k+1)*n, i*n+1:(i+1)*n) = -I;
end
M(k*n+1:(k+1)*n, k*n+1:(k+1)*n) = I;
if strcmp(type, 'whole')
    return
end
M = M(1:k*n, 1:k*n);
end