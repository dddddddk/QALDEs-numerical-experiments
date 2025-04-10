function [U, B] = decompose_backward_matrix(A, d)
% This function generate a decompose of the backward matrix W such that 
% W = U * B
n = length(A); k = length(d);
U = speye(n*k); B = sparse(n*k, n*k);
I = eye(n); 
U(1:n, n+1:2*n) = I;
for i = 3:k
    U(1:n, (i-1)*n+1:i*n) = d(k+3-i)/d(k+2-i) * (-A) * U(1:n, (i-2)*n+1:(i-1)*n) + I;
end
Deno_Matrix = d(k) * (-A) + d(k-1) * I;
for i = (k-2):-1:1
    Deno_Matrix = (-A) * Deno_Matrix + d(i) * I; 
end
B(1:n, (k-1)*n+1:k*n) = Deno_Matrix;
for i = 2:k
    B((i-1)*n+1:i*n, (i-2)*n+1:(i-1)*n) = I;
    B((i-1)*n+1:i*n, (i-1)*n+1:i*n) = d(k+2-i) / d(k+1-i) * A;
end
end