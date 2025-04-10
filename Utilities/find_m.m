function [m, solnorm] = find_m(epsilon, start_m, A, b, x0, T, k, type)
    n = length(A);
    if length(x0) ~= n || length(b) ~= n
        error('Error: The length of x0 or b is wrong!')
    end
    err = 1; m = start_m; solnorm = 0;
    xT = expm(T*A) * x0 + (expm(T*A) - eye(n)) * (A \ b);
    if strcmp(type, 'Taylor')
        while err > epsilon
            m = m + 1;
            h = T / m;
            d = Taylor_Coefficient(k);
            L = construct_m_step_forward_matrix(A*h, d, m);
            vec = construct_vec(x0, h*b, m, d, 'forward');
            z = L \ vec;
            solnorm = norm(z, 2)^2;
            xT_Taylor = z(end-n+1:end);
            err = norm(xT_Taylor - xT, 2) / norm(xT, 2);
        end
        return
    end
    if strcmp(type, 'Pade')
        while err > epsilon
            m = m + 1;
            h = T / m;
            d = Pade_Coefficient(k);
            L = construct_m_step_backward_matrix(A*h, d, m);
            vec = construct_vec(x0, h*b, m, d, 'backward');
            z = L \ vec;
            solnorm = norm(z, 2)^2;
            xT_Pade = z(end-n+1:end);
            err = norm(xT_Pade - xT, 2) / norm(xT, 2);
        end
        return
    end
end