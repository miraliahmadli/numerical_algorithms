function [mu x iter] = power_method(A, x0, epsilon, N)
    x_p = norm(x0, "inf");
    x = x0 / x_p;
    for i = 1:N
        iter = i;
        y = A*x;
        mu = norm(y, "inf");
        
        new_x = y / mu;
        err = norm(x - new_x, "inf") / norm(x, "inf");
        x = new_x;
        if err < epsilon
            break
        end
    end
end


A = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1; 0 3 -1 8];
x0 = [1; 1; 1; 1];
epsilon = 10^(-3);
N = 100;
[mu x iter] = power_method(A, x0, epsilon, N);
disp(sprintf('Power method after %d iterations: mu = %1.4f, x = %1.4f %1.4f %1.4f %1.4f', iter, mu, x));
