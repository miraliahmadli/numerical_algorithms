function [xg Ng] = Randomized_GS(A, b, x0, epsilon, N)
    xg = x0;
    [m, n] = size(A);
    Ng = 0;
    for k = 1:N
        Ng = k;
        r = randperm(n);
        for i = 1:n
            if i == 1
                xg(r(i),1) = 1/A(r(i),r(i)) * (-A(r(i),r(i+1:n))*x0(r(i+1:n),1) + b(r(i)));
            else
                xg(r(i),1) = 1/A(r(i),r(i)) * (-A(r(i),r(1:i-1))*xg(r(1:i-1),1) - A(r(i),r(i+1:n))*x0(r(i+1:n),1) + b(r(i)));
            end
        end
        if norm(xg - x0, inf) / norm(xg, inf) < epsilon
            break
        end
        x0 = xg;
    end
end


A = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1; 0 3 -1 8];
b = [6; 25; -11; 15];
x0 = [0; 0; 0; 0];
epsilon = 10^(-3);
N = 100;
rng(1);
[xg Ng] = Randomized_GS(A, b, x0, epsilon, N);
disp(sprintf('Randomized Gauss-Seidel: x(%d) = %1.4f %1.4f %1.4f %1.4f', Ng, xg));
