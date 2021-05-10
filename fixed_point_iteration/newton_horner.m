% P(x) = x^4 - 4x^2 - 3x + 5
function sol = newton_horner(p0, N, eps)
    function Q = horner(P, pn)
        n = length(P) - 1;
        Q = zeros(n+1, 1);
        Q(1) = P(1);
        for j = 2:n+1
            Q(j) = P(j) + Q(j-1)*pn;
        end
    end
    
    P = [1; 0; -4; -3; 5];
    sol = [p0; 1; 1; 0; -4; -3];
    n = length(P);
    
    for i = 1:N 
        Q = horner(P, p0);
        R = horner(Q, p0);
        
        px0 = Q(5);
        qx0 = R(4);
        sol(1) = p0 - px0/qx0;
        sol(2) = i;
        sol(3) = Q(1);
        sol(4) = Q(2);
        sol(5) = Q(3);
        sol(6) = Q(4);

        if abs(sol(1) - p0) / abs(sol(1)) < eps
            break
        end
        
        p0 = sol(1);
    end
end


p0 = 0.5; N = 10; eps = 0.001;
sol = newton_horner(p0, N, eps);
disp(sprintf('p(%d) = %1.8f', sol(2), sol(1)));
disp(sprintf('Q(x) = %1.2fx^3 + %1.2fx^2 + %1.2fx + %1.2f', sol(3), sol(4), sol(5), sol(6)));
