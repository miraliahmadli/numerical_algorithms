function [ws xs Ns] = SOR(A, b, x0, epsilon, N)
    function [x_next] = next(T, x_old, cj)
        x_next = T*x_old + cj;
    end
    D = diag(diag(A));
    U = (-1) * triu(A, 1);
    L = D - U - A;
    Tj = inv(D)*(D-A); % A = D - L - U
    rho_Tj = max(abs(eig(Tj)));
    ws = 2 / (1+sqrt(1-rho_Tj^2));
    
    Tw = inv(D - ws*L) * ((1-ws)*D + ws*U);

    cw = ws*inv(D-ws*L)*b;
    xs = x0;
    Ns = 0;
    for k = 1:N
        Ns = k;
        xs = next(Tw, x0, cw);
        if norm(xs - x0, inf) / norm(xs, inf) < epsilon
            break
        end
        x0 = xs;
    end
end


A = [10 -1 0 0; -1 11 -1 0; 0 -1 10 -1; 0 0 -1 8];
b = [6; 25; -11; 15];
x0 = [0; 0; 0; 0];
epsilon = 10^(-3);
N = 100;
[ws xs Ns] = SOR(A, b, x0, epsilon, N);
disp(sprintf('SOR(w=%1.4f): x(%d) = %1.4f %1.4f %1.4f %1.4f', ws, Ns, xs));
