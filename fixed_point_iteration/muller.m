% f(x) = x^4 + 2.4x^3 - 12.95x^2 - 34.608x + 91.296 = 0
function sol = muller(p0, p1, p2, N, eps)
    function eval = func(x)
        eval = x^4 + 2.4 * x^3 - 12.95 * x^2 - 34.608 * x + 91.296;
    end
    sol = [p0; 1];
    
    for i = 3:N 
        h1 = p1 - p0;
        h2 = p2 - p1;
        sgm1 = (func(p1) - func(p0)) / h1;
        sgm2 = (func(p2) - func(p1)) / h2;
        d = (sgm2 - sgm1) / (h2 + h1);
        
        b = sgm2 + h2 * d;
        D = (b^2- 4*func(p2)*d)^(1/2);
        if abs(b - D) < abs(b + D)
             E = b + D;
        else
             E = b - D;
        end
        
        h = -2 * func(p2) / E;
        sol(1) = p2 + h;
        sol(2) = i;

        if abs(h) / abs(sol(1)) < eps
            break
        end
        p0 = p1;
        p1 = p2;
        p2 = sol(1);
    end
end


p0 = -1.0; p1 = -2.0; p2 = -3.0; N = 5; eps = 0.001;
sol = muller(p0, p1, p2, N, eps);
disp(sprintf('p(%d) = %1.8f + %1.8fi', sol(2), real(sol(1)), imag(sol(1))));
