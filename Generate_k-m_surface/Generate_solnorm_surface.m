function [error_Matrix, time_Matrix, solnorm_Matrix] = Generate_solnorm_surface(A, x0, b, T, m_list, k_list, type)
n = length(A);
if length(x0) ~= n || length(b) ~= n
    error('Error: The length of x0 or b is wrong!')
end
xT = expm(T*A) * x0 + (expm(T*A) - eye(n)) * (A \ b);
error_Matrix = zeros(length(m_list), length(k_list));
time_Matrix = zeros(length(m_list), length(k_list));
solnorm_Matrix = zeros(length(m_list), length(k_list));
if strcmp(type, 'Taylor')
    for m_index = 1:length(m_list)
        m = m_list(m_index);
        for k_index = 1:length(k_list)
            k = k_list(k_index);
            h = T / m;
            d = Taylor_Coefficient(k);
            L = construct_m_step_forward_matrix(A*h, d, m);
            vec = construct_vec(x0, h*b, m, d, 'forward');
            tic
            z = L \ vec;
            time_Matrix(m_index, k_index) = toc;
            solnorm_Matrix(m_index, k_index) = norm(z, 2)^2;
            xT_Taylor = z(end-n+1:end);
            error_Matrix(m_index, k_index) = log(norm(xT_Taylor - xT, 2) / norm(xT, 2)) / log(10);
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
            vec = construct_vec(x0, h*b, m, d, 'forbackward');
            tic
            z = L \ vec;
            time_Matrix(m_index, k_index) = toc;
            solnorm_Matrix(m_index, k_index) = norm(z, 2)^2;
            xT_Pade = z(end-n+1:end);
            error_Matrix(m_index, k_index) = log(norm(xT_Pade - xT, 2) / norm(xT, 2)) / log(10);
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
            vec = construct_vec(x0, h*b, m, d, 'forbackward');
            tic
            z = L \ vec;
            time_Matrix(m_index, k_index) = toc;
            solnorm_Matrix(m_index, k_index) = norm(z, 2)^2;
            xT_Pade = z(end-n+1:end);
            error_Matrix(m_index, k_index) = log(norm(xT_Pade - xT, 2) / norm(xT, 2)) / log(10);
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
            vec = construct_vec(x0, h*b, m, d, 'backward');
            tic
            z = L \ vec;
            time_Matrix(m_index, k_index) = toc;
            solnorm_Matrix(m_index, k_index) = norm(z, 2)^2;
            xT_Pade = z(end-n+1:end);
            error_Matrix(m_index, k_index) = log(norm(xT_Pade - xT, 2) / norm(xT, 2)) / log(10);
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
            vec = construct_vec(x0, h*b, m, d, 'backward_new');
            tic
            z = L \ vec;
            time_Matrix(m_index, k_index) = toc;
            solnorm_Matrix(m_index, k_index) = norm(z, 2)^2;
            xT_Pade = z(end-n+1:end);
            error_Matrix(m_index, k_index) = log(norm(xT_Pade - xT, 2) / norm(xT, 2)) / log(10);
        end
    end
end
end