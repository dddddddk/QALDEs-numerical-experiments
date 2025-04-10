%% Add Paths
addpath(genpath('./Approximants/'))
addpath(genpath('./Matrix_Construction/'))
addpath(genpath('./Utilities/'))
%% Matrix and Vectors
n = 5;
A = [-2,1,0,0,0;1,-2,1,0,0;0,1,-2,1,0;0,0,1,-2,1;0,0,0,1,-2];
b = ones(n, 1); x0 = ones(n, 1);
T = 1:50;
epsilon = 1e-10;
k_Taylor = 9; k_Pade = 9;
d_Taylor = Taylor_Coefficient(k_Taylor);
d_Pade = Pade_Coefficient(k_Pade);
%% Generate m_Taylor and m_Pade
m_Taylor = zeros(length(T), 1);
m_Pade = zeros(length(T), 1);
cond_Taylor = zeros(length(T), 1);
cond_Pade = zeros(length(T), 1);
succ_Taylor = zeros(length(T), 1);
succ_Pade = zeros(length(T), 1);
solnorm_Taylor = 0; solnorm_Pade = 0;
for i = 1:length(T)
    disp(i)
    xT = expm(T(i)*A) * x0 + (expm(T(i)*A) - eye(n)) * (A \ b);
    xT_norm_square = norm(xT, 2)^2;
    [m_Taylor(i), solnorm_Taylor] = find_m(epsilon, 0, A, b, x0, T(i), k_Taylor, 'Taylor');
    [m_Pade(i), solnorm_Pade] = find_m(epsilon, 0, A, b, x0, T(i), k_Pade, 'Pade');
    succ_Taylor(i) = xT_norm_square / solnorm_Taylor;
    succ_Pade(i) = xT_norm_square / solnorm_Pade;
    h_Taylor = T(i) / m_Taylor(i);
    h_Pade = T(i) / m_Pade(i);
    L_Taylor = construct_m_step_forward_matrix(A*h_Taylor, d_Taylor, m_Taylor(i));
    cond_Taylor(i) = cond(full(L_Taylor));
    clear L_Taylor
    L_Pade = construct_m_step_backward_matrix(A*h_Pade, d_Pade, m_Pade(i));
    cond_Pade(i) = cond(full(L_Pade));
    clear L_Pade
end
%% The T-m plot
figure('Position', [400, 400, 1320, 330])
t = tiledlayout(1, 3);
nexttile
plot(T, [m_Taylor, m_Pade], 'linewidth', 1.1)
set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
legend('Taylor', 'Padé', 'Location', 'northwest')
txt = xlabel('$T$', 'fontsize', 18);
set(txt, 'Interpreter', 'latex');
txt = ylabel('$m^*$', 'fontsize', 18);
set(txt, 'Interpreter', 'latex');
nexttile
plot(T, [cond_Taylor, cond_Pade], 'linewidth', 1.1)
set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
legend('Taylor', 'Padé', 'Location', 'northwest')
txt = xlabel('$T$', 'fontsize', 18);
set(txt, 'Interpreter', 'latex');
txt = ylabel('Condition Number', 'fontsize', 18);
set(txt, 'Interpreter', 'latex');
nexttile
plot(T, [succ_Taylor, succ_Pade], 'linewidth', 1.1)
set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
legend('Taylor', 'Padé', 'Location', 'northeast')
txt = xlabel('$T$', 'fontsize', 18);
set(txt, 'Interpreter', 'latex');
ylabel('Success Probability', 'fontsize', 18)
exportgraphics(t, './fig/fig2.eps', 'ContentType', 'vector')
%% Remove Paths
rmpath(genpath('./Approximants/'))
rmpath(genpath('./Matrix_Construction/'))
rmpath(genpath('./Utilities/'))