%% Add Paths
addpath(genpath('./Approximants/'))
addpath(genpath('./Matrix_Construction/'))
addpath(genpath('./Generate_k-m_surface/'))
%% Set random seed
% rng(2,"twister")
%% Generate matrix and vectors
n = 5; T = 30;
% A special case
A = [-2,1,0,0,0;1,-2,1,0,0;0,1,-2,1,0;0,0,1,-2,1;0,0,0,1,-2];
[U, D] = eig(A);
b = ones(n, 1); x0 = ones(n, 1);
% b = U(:, 1); x0 = U(:, n);
% Diagonalizable matrix with complex eigenvalues with negative real part
% V = randn(n); 
% D = diag(-10 * rand(n, 1) + 1 * 1j * randn(n, 1));
% D = diag(-7 * rand(n, 1) + 7 * 1j * randn(n, 1));
% D = diag(-1 * rand(n, 1) + 10 * 1j * randn(n, 1));
% A = V \ (D * V); 
% Random negative definite matrix
% A = randn(n);
% A = - A' * A;
% Matrix with negative real part eigenvalues
% U = triu(randn(n));
% D = diag(-10 * rand(n, 1) + 1 * 1j * randn(n, 1));
% D = diag(-7 * rand(n, 1) + 7 * 1j * randn(n, 1));
% D = diag(-1 * rand(n, 1) + 10 * 1j * randn(n, 1));
% U = U - diag(diag(U)) + D;
% [Q,~] = qr(randn(n));
% A = Q * U * Q';
% b = randn(n, 1); x0 = randn(n, 1);
fprintf("The eigenvalues of A:\n")
eig_A = eig(A)';
disp(eig_A)
xT = expm(T*A) * x0 + (expm(T*A) - eye(n)) * (A \ b);
xT_norm_square = norm(xT, 2)^2;
m_list = 3:50; k_Taylor = 9; k_Pade = 9;
%% Generate curve about Taylor method
error_Matrix_Taylor = Generate_error_surface(A, x0, b, T, m_list, k_Taylor, 'Taylor');
cond_Matrix_Taylor = Generate_condition_number_surface(A, x0, b, T, m_list, k_Taylor, 'Taylor');
[~,~,solnorm_Matrix_Taylor] = Generate_solnorm_surface(A, x0, b, T, m_list, k_Taylor, 'Taylor');
succ_Taylor = xT_norm_square ./ solnorm_Matrix_Taylor;
%% Generate curve about backward Pade method
[error_Matrix_Pade_back, ~] = Generate_error_surface(A, x0, b, T, m_list, k_Taylor, 'Pade_back');
cond_Matrix_Pade_back = Generate_condition_number_surface(A, x0, b, T, m_list, k_Taylor, 'Pade_back');
[~,~,solnorm_Matrix_Pade_back] = Generate_solnorm_surface(A, x0, b, T, m_list, k_Taylor, 'Pade_back');
succ_Pade = xT_norm_square ./ solnorm_Matrix_Pade_back;
%% Plot
figure('Position', [400, 400, 1320, 330])
% Approximation error
t = tiledlayout(1, 3);
nexttile
plot(m_list, [error_Matrix_Taylor, error_Matrix_Pade_back], 'linewidth', 1.1)
legend('Taylor', 'Padé')
set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
txt = xlabel('$m$', 'fontsize', 18);
set(txt, 'Interpreter', 'latex');
ylabel('Solution Error (log_{10})', 'fontsize', 18)
% Condition number
nexttile
plot(m_list, [log10(cond_Matrix_Taylor), log10(cond_Matrix_Pade_back)], 'linewidth', 1.1)
legend('Taylor', 'Padé')
set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
txt = xlabel('$m$', 'fontsize', 18);
set(txt, 'Interpreter', 'latex');
ylabel('Condtion Number (log_{10})', 'fontsize', 18)
% Success probability
nexttile
plot(m_list, [succ_Taylor, succ_Pade], 'linewidth', 1.1)
set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
legend('Taylor', 'Padé', 'Location', 'northeast')
txt = xlabel('$m$', 'fontsize', 18);
set(txt, 'Interpreter', 'latex');
ylabel('Success Probability', 'fontsize', 18)
exportgraphics(t, './fig/fig1.eps', 'ContentType', 'vector')
%% Remove Paths
rmpath(genpath('./Approximants/'))
rmpath(genpath('./Matrix_Construction/'))
rmpath(genpath('./Generate_k-m_surface/'))