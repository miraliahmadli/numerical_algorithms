function sol = newton(x0, N, epsilon)
    function F = f(x)
        f_1 = 3*x(1) - cos(x(2)*x(3)) - 0.5;
        f_2 = 4*x(1)^2 - 625*x(2)^2 + 2*x(2) - 1;
        f_3 = exp(-x(1)*x(2)) + 20*x(3) + (10*pi - 3)/3;
        F = [f_1; f_2; f_3];
    end
    
    function J = Jacobian(x)
        J = [3, x(3)*sin(x(2)*x(3)), x(2)*sin(x(2)*x(3));
            8*x(1), -1250*x(2) + 2, 0;
            -x(2)*exp(-x(1)*x(2)), -x(1)*exp(-x(1)*x(2)), 20];
    end
    
    sol = [x0; 0];
    xk = x0;
    for k = 1:N
        F = f(xk);
        J = Jacobian(xk);
%         y = -inv(J)*F;
        y = linsolve(J, -F);
        xk = xk + y;
        sol = [xk; k];
        if norm(y, inf) / norm(xk, inf) < epsilon
            break
        end
    end
end
