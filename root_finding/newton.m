% f(x) = sin x - exp(-x)
function sol = newton(p0, N, eps)
    sol = [p0; 1];
    p_old = p0;
    for i = 1:N
        sol(2) = i;
        sol(1) = p_old - ((sin(p_old) - exp(-p_old)) / (cos(p_old) + exp(-p_old)));
        if abs(sol(1) - p_old) / abs(sol(1)) < eps
            break
        end
        p_old = sol(1);
    end
end


p0 = 3.0; N = 5; eps = 0.001;
sol = newton(p0, N, eps);
disp(sprintf('p(%d) = %1.8f', sol(2), sol(1)));
