function [a b c d] = cubic_spline(x, f, opt, g0, gn)
    function [a, b, c, d] = natural_cubic_spline(x, f)
        [n, m] = size(x);
        h = zeros(n-1,1);
        for i = 1:n-1
            h(i) = x(i+1) - x(i);
        end
        
        alpha = zeros(n-2,1);
        for i = 1:n-2
            alpha(i) = 3 / h(i+1) * (f(i+2) - f(i+1)) - 3 / h(i) * (f(i+1) - f(i));
        end
        
        l = zeros(n,1); mu = zeros(n,1); z = zeros(n,1);
        l(1) = 1; l(n) = 1;
        for i = 2:n-1
            l(i) = 2*(x(i+1) - x(i-1)) - h(i-1)*mu(i-1);
            mu(i) = h(i) / l(i);
            z(i) = (alpha(i-1) - h(i-1)*z(i-1)) / l(i);
        end
        
        b = zeros(n-1,1); c = zeros(n,1); d = zeros(n-1,1);
        for j = n-1:-1:1
            c(j) = z(j) - mu(j)*c(j+1);
            b(j) = (f(j+1) - f(j)) / h(j) - h(j)*(c(j+1) + 2*c(j))/3;
            d(j) = (c(j+1) - c(j)) / (3*h(j));
        end
        c(end) = [];
        a = f(1:n-1);
    end
    
    function [a, b, c, d] = clamped_cubic_spline(x, f, g0, gn)
        [n, m] = size(x);
        h = zeros(n-1,1);
        for i = 1:n-1
            h(i) = x(i+1) - x(i);
        end
        
        alpha = zeros(n,1);
        alpha(1) = 3*(f(2) - f(1))/ h(1) - 3*g0;
        alpha(n) = 3*gn - 3*(f(n) - f(n-1))/ h(n-1);
        for i = 2:n-1
            alpha(i) = 3 / h(i) * (f(i+1) - f(i)) - 3 / h(i-1) * (f(i) - f(i-1));
        end
        
        l = zeros(n,1); mu = zeros(n,1); z = zeros(n,1);
        l(1) = 2*h(1); mu(1) = 0.5; z(1) = alpha(1) / l(1);
        for i = 2:n-1
            l(i) = 2*(x(i+1) - x(i-1)) - h(i-1)*mu(i-1);
            mu(i) = h(i) / l(i);
            z(i) = (alpha(i) - h(i-1)*z(i-1)) / l(i);
        end
        l(n) = h(n-1)*(2-mu(n-1));
        z(n) = (alpha(n) - h(n-1)*z(n-1))/l(n);
        
        b = zeros(n-1,1); c = zeros(n,1); d = zeros(n-1,1);
        c(n) = z(n);
        for j = n-1:-1:1
            c(j) = z(j) - mu(j)*c(j+1);
            b(j) = (f(j+1) - f(j)) / h(j) - h(j)*(c(j+1) + 2*c(j))/3;
            d(j) = (c(j+1) - c(j)) / (3*h(j));
        end
        c(end) = [];
        a = f(1:n-1);
    end
    
    if opt == 0
        [a, b, c, d] = natural_cubic_spline(x, f);
    else
        [a, b, c, d] = clamped_cubic_spline(x, f, g0, gn);
    end
end


x = [0; 1; 2; 3];
f = [1; exp(1); exp(2); exp(3)];
disp('===Natural cubic spline===');
[a b c d] = cubic_spline(x, f, 0);
disp(sprintf('a = [%1.5f %1.5f %1.5f]', a));
disp(sprintf('b = [%1.5f %1.5f %1.5f]', b));
disp(sprintf('c = [%1.5f %1.5f %1.5f]', c));
disp(sprintf('d = [%1.5f %1.5f %1.5f]', d));
