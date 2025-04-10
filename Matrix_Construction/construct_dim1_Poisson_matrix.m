function T = construct_dim1_Poisson_matrix(n, bdycond)
% This function construct a finite difference matrix corresponding to
% d^2/dx^2
T = sparse(n, n);
if strcmp(bdycond, 'Dirichlet')
    T(1, 1) = -2;
    for i = 2:n
        T(i, i) = -2;
        T(i-1, i) = 1; T(i, i-1) = 1;
    end
end
end