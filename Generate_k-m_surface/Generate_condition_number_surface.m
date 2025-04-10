function cond_Matrix = Generate_condition_number_surface(A, x0, b, T, m_list, k_list, type)
n = length(A);
if length(x0) ~= n || length(b) ~= n
    error('Error: The length of x0 or b is wrong!')
end
cond_Matrix = zeros(length(m_list), length(k_list));
if strcmp(type, 'Taylor')
    for m_index = 1:length(m_list)
        m = m_list(m_index);
        for k_index = 1:length(k_list)
            k = k_list(k_index);
            h = T / m;
            d = Taylor_Coefficient(k);
            L = construct_m_step_forward_matrix(A*h, d, m);
            cond_Matrix(m_index, k_index) = cond(full(L));
        end
    end
end
if strcmp(type, 'Pade')
    for m_index = 1:length(m_list)
        m = m_list(m_index);
        for k_index = 1:length(k_list)
            k = k_list(k_index);
            h = T / m;
            d = Pade_Coefficient(k);
            L = construct_m_step_forbackward_matrix(A*h, d, m);
            cond_Matrix(m_index, k_index) = cond(full(L));
        end
    end
end
if strcmp(type, 'Pade_new')
    for m_index = 1:length(m_list)
        m = m_list(m_index);
        for k_index = 1:length(k_list)
            k = k_list(k_index);
            h = T / m;
            d = Pade_Coefficient(k);
            L = construct_m_step_forbackward_matrix_new(A*h, d, m);
            cond_Matrix(m_index, k_index) = cond(full(L));
        end
    end
end
if strcmp(type, 'Pade_back')
    for m_index = 1:length(m_list)
        m = m_list(m_index);
        for k_index = 1:length(k_list)
            k = k_list(k_index);
            h = T / m;
            d = Pade_Coefficient(k);
            L = construct_m_step_backward_matrix(A*h, d, m);
            cond_Matrix(m_index, k_index) = cond(full(L));
        end
    end
end
if strcmp(type, 'Pade_back_new')
    for m_index = 1:length(m_list)
        m = m_list(m_index);
        for k_index = 1:length(k_list)
            k = k_list(k_index);
            h = T / m;
            d = Pade_Coefficient(k);
            L = construct_m_step_backward_matrix_new(A*h, d, m);
            cond_Matrix(m_index, k_index) = cond(full(L));
        end
    end
end
if strcmp(type, 'Pade_deno')
    for m_index = 1:length(m_list)
        m = m_list(m_index);
        for k_index = 1:length(k_list)
            k = k_list(k_index);
            h = T / m;
            d = Pade_Coefficient(k);
            L = construct_Denominator(A*h, d);
            cond_Matrix(m_index, k_index) = cond(L);
        end
    end
end
end