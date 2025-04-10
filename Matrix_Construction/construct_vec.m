function vec = construct_vec(x, b, m, d, type)
% If type is 'forward', then this function construct the vector
% corresponding to the m step forward matrix. If type is 'forbackward',
% then this function return the vector corresponding to the m step
% forward-backward matrix. If type is 'backward', then this function return
% the vector corresponding to the m step backward matrix. 
% x, b: given vectors
% m: time steps
% d: the coefficient list
n = length(x); k = length(d);
if length(b) ~= n
    error('Error: The size of x0 and b is different!')
end
if strcmp(type, 'forward')
%     vec = zeros((1+k*m)*n, 1);
    vec = sparse((1+k*m)*n, 1);
    vec(1:n) = x;
    vec(n+1:2*n) = d(2) * b;
    sub_vec = sparse(k*n, 1); sub_vec(n+1:2*n) = d(2) * b;
    for i = 2:m
        vec((i-1)*n*k+1:i*k*n) = sub_vec;
    end
end
if strcmp(type, 'forbackward')
    vec = sparse((m+1)*k*n, 1);
%     vec = zeros((m+1)*k*n, 1);
    vec(1:n) = x; vec(n+1:2*n) = d(2) * b;
    sub_vec = sparse(k*n, 1); sub_vec((k-1)*n+1:k*n) = -d(2) * b;
    for i = 1:m
        vec(i*k*n+1:(i+1)*n*k) = sub_vec;
    end
end
if strcmp(type, 'backward')
    vec = sparse((m*k + 1)*n, 1);
    vec(1:n) = x; vec((k-1)*n+1:k*n) = -d(2) * b;
    sub_vec = sparse(k*n, 1); sub_vec((k-1)*n+1:k*n) = -d(2) * b;
    for i = 1:m-1
        vec(i*k*n+1:(i+1)*n*k) = sub_vec;
    end
end
if strcmp(type, 'backward_new')
    vec = sparse((m*k + 1)*n, 1);
    vec(1:n) = x / sqrt(k); vec((k-1)*n+1:k*n) = -d(2) * b;
    sub_vec = sparse(k*n, 1); sub_vec((k-1)*n+1:k*n) = -d(2) * b;
    for i = 1:m-1
        vec(i*k*n+1:(i+1)*n*k) = sub_vec;
    end
end
end