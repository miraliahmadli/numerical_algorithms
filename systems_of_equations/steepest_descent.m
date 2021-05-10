function [xs Ns] = steepest_descent(A, b, x0, epsilon, N)
    function [x_next] = next(x0, A, r)
        tk = dot(r,r) / dot(r, A*r);
        x_next = x0 + tk * r;
    end

    r0 = b - A*x0;
    xs = x0;
    Ns = 0;
    for k = 1:N
        Ns = k;
        xs = next(x0, A, r0);
        if norm(xs - x0, inf) / norm(xs, inf) < epsilon
            break
        end
        x0 = xs;
        r0 = b - A*x0;
    end
end


A = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1; 0 3 -1 8];
b = [6; 25; -11; 15];
x0 = [0; 0; 0; 0];
epsilon = 10^(-3);
N = 100;
[xs Ns] = steepest_descent(A, b, x0, epsilon, N);
disp(sprintf('Steepest descent: x(%d) = %1.4f %1.4f %1.4f %1.4f', Ns, xs));
