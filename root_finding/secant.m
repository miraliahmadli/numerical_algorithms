% f(x) = sin x - exp(-x)
function sol = secant(p0, p1, N, eps)
    sol = [p1; 1];
    p_old = p1;
    p_old2 = p0;
    for i = 1:N
        sol(2) = i + 1;
        sol(1) = p_old - (((sin(p_old) - exp(-p_old)) * (p_old - p_old2)) / ((sin(p_old) - exp(-p_old)) - (sin(p_old2) - exp(-p_old2))));
        if abs(sol(1) - p_old) / abs(sol(1)) < eps
            break
        end
        p_old2 = p_old;
        p_old = sol(1);
    end
end


p0 = 3.0; p1 = 3.5; N = 5; eps = 0.001;
sol = secant(p0, p1, N, eps);
disp(sprintf('p(%d) = %1.8f', sol(2), sol(1)));
