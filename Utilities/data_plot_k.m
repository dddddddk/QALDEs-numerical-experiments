function data_plot_k(N, k_Taylor, k_Pade, cond_Taylor, cond_Pade, solnorm_Taylor, solnorm_Pade, save_file)
    k_Taylor = k_Taylor(:, 1:N);
    k_Pade = k_Pade(:, 1:N);
    cond_Taylor = cond_Taylor(:, 1:N);
    cond_Pade = cond_Pade(:, 1:N);
    solnorm_Taylor = solnorm_Taylor(:, 1:N);
    solnorm_Pade = solnorm_Pade(:, 1:N);
    eps_list = 1:14;
    %% Compute sample mean and standard deviation
    % Sample mean
    mean_k_Taylor = k_Taylor * ones(N, 1) / N;
    mean_k_Pade = k_Pade * ones(N, 1) / N;
    mean_cond_Taylor = cond_Taylor * ones(N, 1) / N;
    mean_cond_Pade = cond_Pade * ones(N, 1) / N;
    mean_solnorm_Taylor = solnorm_Taylor * ones(N, 1) / N;
    mean_solnorm_Pade = solnorm_Pade * ones(N, 1) / N;
    % Standard deviation
    std_k_Taylor = std(k_Taylor, 0, 2);
    std_k_Pade = std(k_Pade, 0, 2);
    std_cond_Taylor = std(cond_Taylor, 0, 2);
    std_cond_Pade = std(cond_Pade, 0, 2);
    std_solnorm_Taylor = std(solnorm_Taylor, 0, 2);
    std_solnorm_Pade = std(solnorm_Pade, 0, 2);
    %% Plot
    figure('Position', [400, 400, 880, 330])
%     figure('Position', [400, 400, 1200, 300])
    % The eps-k plot
    t = tiledlayout(1, 2);
%     subplot(1, 3, 1)
    nexttile
    hold on
    fill([eps_list, eps_list(end:-1:1)], [mean_k_Taylor'- std_k_Taylor', mean_k_Taylor(end:-1:1)'+ std_k_Taylor(end:-1:1)'], 'r', 'FaceColor', '#0072BD', 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Area1');
    fill([eps_list, eps_list(end:-1:1)], [mean_k_Pade'- std_k_Pade', mean_k_Pade(end:-1:1)'+ std_k_Pade(end:-1:1)'], 'r', 'FaceColor', '#E27E82', 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Area2');
    h = plot(eps_list, [mean_k_Taylor, mean_k_Pade], 'linewidth', 1.1);
    xticks([1 5 10 14])
    xticklabels({'10^{-1}', '10^{-5}', '10^{-10}', '10^{-14}'})
    set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times');
    legend(h, 'Taylor', 'Padé', 'Location', 'northwest')
    xlabel('Relative Solution Error', 'fontsize', 18)
    txt = ylabel('$k^*$', 'fontsize', 18);
    set(txt, 'Interpreter', 'latex');
    % The eps-condition_number plot
%     subplot(1, 3, 2)
    nexttile
    hold on
    fill([eps_list, eps_list(end:-1:1)], [mean_cond_Taylor'- std_cond_Taylor', mean_cond_Taylor(end:-1:1)'+ std_cond_Taylor(end:-1:1)'], 'r', 'FaceColor', '#0072BD', 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Area3');
    fill([eps_list, eps_list(end:-1:1)], [mean_cond_Pade'- std_cond_Pade', mean_cond_Pade(end:-1:1)'+ std_cond_Pade(end:-1:1)'], 'r', 'FaceColor', '#E27E82', 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Area4');
    h = plot(eps_list, [mean_cond_Taylor, mean_cond_Pade], 'linewidth', 1.1);
    xticks([1 5 10 14])
    xticklabels({'10^{-1}', '10^{-5}', '10^{-10}', '10^{-14}'})
    set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
    legend(h, 'Taylor', 'Padé', 'Location', 'northwest')
    xlabel('Relative Solution Error', 'fontsize', 18)
    ylabel('Condition Number', 'fontsize', 18)
    % The eps-solution norm plot
%     subplot(1, 3, 3)
%     nexttile
%     h1 = plot(eps_list, [mean_solnorm_Taylor, mean_solnorm_Pade], 'linewidth', 1.1);
%     set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
%     hold on
%     h2 = fill([eps_list, eps_list(end:-1:1)], [mean_solnorm_Taylor'- std_solnorm_Taylor', mean_solnorm_Taylor(end:-1:1)'+ std_solnorm_Taylor(end:-1:1)'], 'r');
%     set(h2, 'FaceColor', '#0072BD', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
%     h3 = fill([eps_list, eps_list(end:-1:1)], [mean_solnorm_Pade'- std_solnorm_Pade', mean_solnorm_Pade(end:-1:1)'+ std_solnorm_Pade(end:-1:1)'], 'r');
%     set(h3, 'FaceColor', '#E27E82', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
%     hold off
%     legend(h1, 'Taylor', 'Padé', 'Location', 'northwest')
%     xlabel('Relative Solution Error (-log_{10})')
%     ylabel('Solution Norm Suqare')
    exportgraphics(t, save_file, 'ContentType','vector')
end