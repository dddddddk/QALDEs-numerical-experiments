function T = construct_T(n, k, type)
% n: the block size of T
% k: the number of blocks in T
% type: The type of T matrix, '1' or '2'
I = eye(n);
T = sparse(k*n, k*n);
if strcmp(type, '1')
    for i = 1:k
        T(1:n, (i-1)*n+1:i*n) = -I;
    end
end
if strcmp(type, '2') 
    for i = 1:k
        T(1:n, (i-1)*n+1:i*n) = (-1)^(k+1-i) * I;
    end
end
end