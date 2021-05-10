function [xg Ng] = GS(A, b, x0, epsilon, N)
    xg = x0;
    [m, n] = size(A);
    Ng = 0;
    for k = 1:N
        Ng = k;
        for i = 1:n
            if i == 1
                xg(i,1) = 1/A(i,i) * (-A(i,i+1:n)*x0(i+1:n,1) + b(i));
            else
                xg(i,1) = 1/A(i,i) * (-A(i,1:i-1)*xg(1:i-1,1) - A(i,i+1:n)*x0(i+1:n,1) + b(i));
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
[xg Ng] = GS(A, b, x0, epsilon, N);
disp(sprintf('Gauss-Seidel: x(%d) = %1.4f %1.4f %1.4f %1.4f', Ng, xg));
