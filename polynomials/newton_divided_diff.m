function [a b] = newton_divided_diff(x, f)
    function coeffs = newton_method(x, f)
        [n, m] = size(x);
        F = zeros(n,n);
        F(:,1) = f;
        for i = 1:n-1
            for j = 1:i
                F(i+1,j+1) = (F(i+1,j) - F(i, j)) / (x(i+1) - x(i-j+1));
            end
        end
        coeffs = diag(F);
    end

    a = newton_method(x, f);
    
    x = flip(x);
    f = flip(f);
    b = newton_method(x, f);
end


x = [0.0; 0.2; 0.4; 0.6; 0.8];
f = [1.0000; 1.22140; 1.49182; 1.82212; 2.22554];
[a b] = newton_divided_diff(x, f);
disp(sprintf('a = [%1.5f %1.5f %1.5f %1.5f %1.5f]', a));
disp(sprintf('b = [%1.5f %1.5f %1.5f %1.5f %1.5f]', b));
