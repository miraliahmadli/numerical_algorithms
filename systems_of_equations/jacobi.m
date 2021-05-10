function [xj Nj] = Jacobi(A, b, x0, epsilon, N)
    function [x_next] = next(T, x_old, cj)
        x_next = T*x_old + cj;
    end
    D = diag(diag(A));
    T = inv(D)*(D-A); % A = D - L - U
    cj = inv(D)*b;
    xj = x0;
    Nj = 0;
    for k = 1:N
        Nj = k;
        xj = next(T, x0, cj);
        if norm(xj - x0, inf) / norm(xj, inf) < epsilon
            break
        end
        x0 = xj;
    end
end
