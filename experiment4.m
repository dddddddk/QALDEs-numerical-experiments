%% Add Paths
addpath(genpath('./Approximants/'))
addpath(genpath('./Matrix_Construction/'))
addpath(genpath('./Utilities/'))
%% Matrix and Vectors
n = 5; N = 100;
b = ones(n, 1); x0 = ones(n, 1);
T = 1; m = 1;
eps_list = 10.^(-(1:14));
%% Generate m_Taylor and m_Pade
k_Taylor = zeros(length(eps_list), N);
k_Pade = zeros(length(eps_list), N);
cond_Taylor = zeros(length(eps_list), N);
cond_Pade = zeros(length(eps_list), N);
solnorm_Taylor = zeros(length(eps_list), N);
solnorm_Pade = zeros(length(eps_list), N);
for j = 1:N
    % Matrix with negative real part eigenvalues
    U = triu(randn(n));
    D = 5 * diag(-abs(randn(n, 1)) + 1j * randn(n, 1));
    U = U - diag(diag(U)) + D;
    [Q,~] = qr(randn(n));
    A = Q * U * Q';
    A = A / norm(A, 2);
    xT = expm(T*A) * x0 + (expm(T*A) - eye(n)) * (A \ b);
    for i = 1:length(eps_list)
        epsilon = eps_list(i);
        k_Taylor(i, j) = find_k(epsilon, A, b, x0, T, m, 'Taylor');
        d_Taylor = Taylor_Coefficient(k_Taylor(i, j));
        k_Pade(i, j) = find_k(epsilon, A, b, x0, T, m, 'Pade');
        d_Pade = Pade_Coefficient(k_Pade(i, j));
        L_Taylor = construct_m_step_forward_matrix(A, d_Taylor, m);
        vec_Taylor = construct_vec(x0, b, m, d_Taylor, 'forward');
        z_Taylor = L_Taylor \ vec_Taylor;
        solnorm_Taylor(i, j) = norm(z_Taylor, 2)^2;
        xT_Taylor = z_Taylor(end-n+1:end);
        cond_Taylor(i, j) = cond(full(L_Taylor));
        clear L_Taylor
        L_Pade = construct_m_step_backward_matrix(A, d_Pade, m);
        vec_Pade = construct_vec(x0, b, m, d_Pade, 'backward');
        z_Pade = L_Pade \ vec_Pade;
        solnorm_Pade(i, j) = norm(z_Pade, 2)^2;
        xT_Pade = z_Pade(end-n+1:end);
        cond_Pade(i, j) = cond(full(L_Pade));
        clear L_Pade
    end
    fprintf(['j:', num2str(j), '\n']);
end
data_plot_k(N, k_Taylor, k_Pade, cond_Taylor, cond_Pade, solnorm_Taylor, solnorm_Pade, './fig/fig4.eps');
drawnow;
%% Remove Paths
rmpath(genpath('./Approximants/'))
rmpath(genpath('./Matrix_Construction/'))
rmpath(genpath('./Utilities/'))