%% Add Paths
addpath(genpath('./Approximants/'))
addpath(genpath('./Matrix_Construction/'))
addpath(genpath('./Utilities/'))
%% Matrix and Vectors
n = 5; N = 100;
b = ones(n, 1); x0 = ones(n, 1);
T = 1:50;
epsilon = 1e-10;
k_Taylor = 9; k_Pade = 9;
d_Taylor = Taylor_Coefficient(k_Taylor); d_Pade = Pade_Coefficient(k_Pade);
%% Generate m_Taylor and m_Pade
m_Taylor = zeros(length(T), N);
m_Pade = zeros(length(T), N);
cond_Taylor = zeros(length(T), N);
cond_Pade = zeros(length(T), N);
succ_Taylor = zeros(length(T), N);
succ_Pade = zeros(length(T), N);
solnorm_Taylor = 0; solnorm_Pade = 0;
for j = 1:N
    % Matrix with negative real part eigenvalues
    U = triu(randn(n));
    D = diag(-abs(randn(n, 1)) + 10 * 1j * randn(n, 1));
    U = U - diag(diag(U)) + D;
    [Q,~] = qr(randn(n));
    A = Q * D * Q';
    A = A / norm(A, 2);
    fprintf('i: ')
    for i = 1:length(T)
        fprintf([num2str(i), ' ']);
        if i == 1
            m_start_Taylor = 0;
            m_start_Pade = 0;
        else
            m_start_Taylor = m_Taylor(i-1, j) - 1;
            m_start_Pade = m_Pade(i-1, j) - 1;
        end
        xT = expm(T(i)*A) * x0 + (expm(T(i)*A) - eye(n)) * (A \ b);
        xT_norm_square = norm(xT, 2)^2;
        [m_Taylor(i, j), solnorm_Taylor] = find_m(epsilon, m_start_Taylor, A, b, x0, T(i), k_Taylor, 'Taylor');
        [m_Pade(i, j), solnorm_Pade] = find_m(epsilon, m_start_Pade, A, b, x0, T(i), k_Pade, 'Pade');
        succ_Taylor(i,j) = xT_norm_square / solnorm_Taylor;
        succ_Pade(i,j) = xT_norm_square / solnorm_Pade;
        h_Taylor = T(i) / m_Taylor(i, j);
        h_Pade = T(i) / m_Pade(i, j);
        L_Taylor = construct_m_step_forward_matrix(A*h_Taylor, d_Taylor, m_Taylor(i, j));
        cond_Taylor(i, j) = cond(full(L_Taylor));
        clear L_Taylor
        L_Pade = construct_m_step_backward_matrix(A*h_Pade, d_Pade, m_Pade(i, j));
        cond_Pade(i, j) = cond(full(L_Pade));
        clear L_Pade
    end
    fprintf('\n')
    fprintf(['j:', num2str(j), '\n']);
    if mod(j, 1) == 0
        % Write in .csv
        writematrix(m_Taylor, './data/exp3/m_Taylor.csv');
        writematrix(m_Pade, './data/exp3/m_Pade.csv');
        writematrix(cond_Taylor, './data/exp3/cond_Taylor.csv');
        writematrix(cond_Pade, './data/exp3/cond_Pade.csv');
        writematrix(succ_Taylor, './data/exp3/succ_Taylor.csv');
        writematrix(succ_Pade, './data/exp3/succ_Pade.csv');
        %% Plot
        close
        data_plot(j, m_Taylor, m_Pade, cond_Taylor, cond_Pade, succ_Taylor, succ_Pade, './fig/fig3.eps');
        drawnow;
    end
end
%% Remove Paths
rmpath(genpath('./Approximants/'))
rmpath(genpath('./Matrix_Construction/'))
rmpath(genpath('./Utilities/'))