function [xp Np] = preconditioned_conjugate_gradient(A, b, x0, epsilon, N)
    C = diag(sqrt(diag(A)));
    C_inv = inv(C);
    r = b - A*x0; 
    w = C_inv*r;
    v = C_inv*w;
    alpha = dot(w,w);

    xp = x0;
    Np = 1;
    for k = 1:N
        Np = k+1;
        if norm(v, inf) < epsilon
            break
        end
        
        u = A*v;
        t = alpha / dot(v,u);
        xp = xp + t*v;
        r = r - t*u;
        w = C_inv*r;
        beta = dot(w,w);
        
        if abs(beta) < epsilon
            if norm(r, inf) < epsilon
                break
            end
        end

        s = beta / alpha;
        v = C_inv*w + s*v;
        alpha = beta;
    end
end


A = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1; 0 3 -1 8];
b = [6; 25; -11; 15];
x0 = [0; 0; 0; 0];
epsilon = 10^(-3);
N = 100;
[xp Np] = preconditioned_conjugate_gradient(A, b, x0, epsilon, N);
disp(sprintf('Preconditioned conjugate gradient: x(%d) = %1.4f %1.4f %1.4f %1.4f', Np, xp));
