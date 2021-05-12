function [a b] = divided_diff_for_Hermite(x, f, g)
    function coeffs = Hermite(x, f, g)
        [n, m] = size(x);
        z = zeros(2*n);
        Q = zeros(2*n,2*n);
        for i = 0:n-1
            z(2*i + 1) = x(i + 1);
            z(2*i + 2) = x(i + 1);
            Q(2*i+1, 1) = f(i+1);
            Q(2*i+2, 1) = f(i+1);
            Q(2*i+2, 2) = g(i+1);
            if i ~= 0
                Q(2*i+1, 2) = (Q(2*i+1, 1) - Q(2*i, 1)) / (z(2*i+1) - z(2*i));
            end
        end
        
        for i = 2:2*n-1
            for j = 2:i
                Q(i+1, j+1) = (Q(i+1,j) - Q(i, j)) / (z(i+1) - z(i-j+1));
            end
        end
        coeffs = diag(Q);
    end

    a = Hermite(x, f, g);
    
    x = flip(x);
    f = flip(f);
    g = flip(g);
    b = Hermite(x, f, g);
end


x = [8.0; 8.6];
f = [17.56492; 18.50515];
g = [3.116256; 3.151762];
[a b] = divided_diff_for_Hermite(x, f, g);
disp(sprintf('a = [%1.5f %1.5f %1.5f %1.5f]', a));
disp(sprintf('b = [%1.5f %1.5f %1.5f %1.5f]', b));
