function data_plot(N, m_Taylor, m_Pade, cond_Taylor, cond_Pade, succ_Taylor, succ_Pade, save_file)
    m_Taylor = m_Taylor(:, 1:N);
    m_Pade = m_Pade(:, 1:N);
    cond_Taylor = cond_Taylor(:, 1:N);
    cond_Pade = cond_Pade(:, 1:N);
    succ_Taylor = succ_Taylor(:, 1:N);
    succ_Pade = succ_Pade(:, 1:N);
    T = 1:50;
    %% Compute sample mean and standard deviation
    % Sample mean
    mean_m_Taylor = m_Taylor * ones(N, 1) / N;
    mean_m_Pade = m_Pade * ones(N, 1) / N;
    mean_cond_Taylor = cond_Taylor * ones(N, 1) / N;
    mean_cond_Pade = cond_Pade * ones(N, 1) / N;
    mean_succ_Taylor = succ_Taylor * ones(N, 1) / N;
    mean_succ_Pade = succ_Pade * ones(N, 1) / N;
    % Standard deviation
    std_m_Taylor = std(m_Taylor, 0, 2);
    std_m_Pade = std(m_Pade, 0, 2);
    std_cond_Taylor = std(cond_Taylor, 0, 2);
    std_cond_Pade = std(cond_Pade, 0, 2);
    std_succ_Taylor = std(succ_Taylor, 0, 2);
    std_succ_Pade = std(succ_Pade, 0, 2);
    %% Plot
    figure('Position', [400, 400, 1320, 330])
    % The T-m plot
    t = tiledlayout(1, 3);
    nexttile
    hold on
    h2 = fill([T, T(end:-1:1)], [mean_m_Taylor'- std_m_Taylor', mean_m_Taylor(end:-1:1)'+ std_m_Taylor(end:-1:1)'], 'r');
    set(h2, 'FaceColor', '#0072BD', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    h3 = fill([T, T(end:-1:1)], [mean_m_Pade'- std_m_Pade', mean_m_Pade(end:-1:1)'+ std_m_Pade(end:-1:1)'], 'r');
    set(h3, 'FaceColor', '#E27E82', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    h1 = plot(T, [mean_m_Taylor, mean_m_Pade], 'linewidth', 1.1);
    set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
    hold off
    legend(h1, 'Taylor', 'Padé', 'Location', 'northwest')
    txt = xlabel('$T$', 'fontsize', 18);
    set(txt, 'Interpreter', 'latex');
    txt = ylabel('$m^*$', 'fontsize', 18);
    set(txt, 'Interpreter', 'latex');
    % The T-condition_number plot
    nexttile
    hold on
    h2 = fill([T, T(end:-1:1)], [mean_cond_Taylor'- std_cond_Taylor', mean_cond_Taylor(end:-1:1)'+ std_cond_Taylor(end:-1:1)'], 'r');
    set(h2, 'FaceColor', '#0072BD', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    h3 = fill([T, T(end:-1:1)], [mean_cond_Pade'- std_cond_Pade', mean_cond_Pade(end:-1:1)'+ std_cond_Pade(end:-1:1)'], 'r');
    set(h3, 'FaceColor', '#E27E82', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    h1 = plot(T, [mean_cond_Taylor, mean_cond_Pade], 'linewidth', 1.1);
    set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
    hold off
    legend(h1, 'Taylor', 'Padé', 'Location', 'northwest')
    txt = xlabel('$T$', 'fontsize', 18);
    set(txt, 'Interpreter', 'latex');
    ylabel('Condition Number', 'fontsize', 18);
    % The T-success probability plot
%     subplot(1, 3, 3)
    nexttile
    hold on
    h2 = fill([T, T(end:-1:1)], [mean_succ_Taylor'- std_succ_Taylor', mean_succ_Taylor(end:-1:1)'+ std_succ_Taylor(end:-1:1)'], 'r');
    set(h2, 'FaceColor', '#0072BD', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    h3 = fill([T, T(end:-1:1)], [mean_succ_Pade'- std_succ_Pade', mean_succ_Pade(end:-1:1)'+ std_succ_Pade(end:-1:1)'], 'r');
    set(h3, 'FaceColor', '#E27E82', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    h1 = plot(T, [mean_succ_Taylor, mean_succ_Pade], 'linewidth', 1.1);
    set(gca, 'linewidth', 1.1, 'fontsize', 12, 'fontname', 'times')
    hold off
    legend(h1, 'Taylor', 'Padé', 'Location', 'northeast')
    txt = xlabel('$T$', 'fontsize', 18);
    set(txt, 'Interpreter', 'latex');
    ylabel('Success Probability', 'fontsize', 18);
    exportgraphics(t, save_file, 'ContentType', 'vector')
end