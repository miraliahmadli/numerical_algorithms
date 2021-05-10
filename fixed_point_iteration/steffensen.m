% x = g(x) = sqrt(exp(x) / 3)
function sol = steffensen(p0, N, eps)
    sol = [p0; 1];
    for i = 1:N
        p1 = sqrt(exp(p0) / 3);
        p2 = sqrt(exp(p1) / 3);
        sol(2) = i;
        sol(1) = p0 - ((p1 - p0)^2 / (p2 - 2*p1 + p0));
        if abs(sol(1) - p0) / abs(sol(1)) < eps
            break
        end
        p0 = sol(1);
    end
end


p0 = -2.0; N = 5; eps = 0.001;
sol = steffensen(p0, N, eps);
disp(sprintf('p(%d) = %1.8f', sol(2), sol(1)));
