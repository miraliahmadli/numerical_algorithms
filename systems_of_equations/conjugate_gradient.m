function [xc Nc] = conjugate_gradient(A, b, x0, epsilon, N)
    function [x_next, r_next, v_next] = next(x0, A, r, v)
        tk = dot(r,r) / dot(v, A*v);
        x_next = x0 + tk * v;
        
        r_next = r - tk * A * v;
        
        sk = dot(r_next, r_next) / dot(r, r);
        v_next = r_next + sk*v;
    end

    r0 = b - A*x0; 
    v1 = r0;
    xc = x0;
    Nc = 0;
    for k = 1:N
        Nc = k;
        [xc, rc, vc] = next(x0, A, r0, v1);
        if norm(xc - x0, inf) / norm(xc, inf) < epsilon
            break
        end
        x0 = xc;
        r0 = rc;
        v1 = vc;
    end
end


A = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1; 0 3 -1 8];
b = [6; 25; -11; 15];
x0 = [0; 0; 0; 0];
epsilon = 10^(-3);
N = 100;
[xc Nc] = conjugate_gradient(A, b, x0, epsilon, N);
disp(sprintf('Conjugate gradient: x(%d) = %1.4f %1.4f %1.4f %1.4f', Nc, xc));
