function vec = construct_block_vec(X, B, m, d, type)
% If type is 'forward', then this function construct the vector
% corresponding to the m step forward matrix. If type is 'forbackward',
% then this function return the vector corresponding to the m step
% forward-backward matrix. If type is 'backward', then this function return
% the vector corresponding to the m step backward matrix. 
% x, b: given vectors
% m: time steps
% d: the coefficient list
size_X = size(X); size_B = size(B);
n_row = size_X(1); n_col = size_X(2);
k = length(d);
if size_B(1) ~= n_row
    error('Error: The size of x0 and b is different!')
end
if size_B(2) ~= n_col
    error('Error: The size of x0 and b is different!')
end
if strcmp(type, 'forward')
    vec = sparse((1+k*m)*n_row, n_col);
    vec(1:n_row, :) = X;
    vec(n_row+1:2*n_row, :) = d(2) * B;
    sub_vec = sparse(k*n_row, n_col); sub_vec(n_row+1:2*n_row, :) = d(2) * B;
    for i = 2:m
        vec((i-1)*n_row*k+1:i*k*n_row, :) = sub_vec;
    end
end
if strcmp(type, 'forbackward')
    vec = sparse((m+1)*k*n_row, n_col);
    vec(1:n_row, :) = X; vec(n_row+1:2*n_row, :) = d(2) * B;
    sub_vec = sparse(k*n_row, n_col); sub_vec((k-1)*n_row+1:k*n_row) = -d(2) * B;
    for i = 1:m
        vec(i*k*n_row+1:(i+1)*n_row*k, :) = sub_vec;
    end
end
if strcmp(type, 'backward')
    vec = sparse((m*k + 1)*n_row, n_col);
    vec(1:n_row, :) = X; vec((k-1)*n_row+1:k*n_row, :) = -d(2) * B;
    sub_vec = sparse(k*n_row, n_col); sub_vec((k-1)*n_row+1:k*n_row, :) = -d(2) * B;
    for i = 1:m-1
        vec(i*k*n_row+1:(i+1)*n_row*k, :) = sub_vec;
    end
end
if strcmp(type, 'backward_new')
    vec = sparse((m*k + 1)*n_row, n_col);
    vec(1:n_row, :) = X / sqrt(k); vec((k-1)*n_row+1:k*n_row, :) = -d(2) * B;
    sub_vec = sparse(k*n_row, n_col); sub_vec((k-1)*n_row+1:k*n_row, :) = -d(2) * B;
    for i = 1:m-1
        vec(i*k*n_row+1:(i+1)*n_row*k, :) = sub_vec;
    end
end
end