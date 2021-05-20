function [s1 s2 s3] = composite_integration(f, a, b, n)
    function s = simpson(f, a, b, n)
        h = (b-a)/n;
        x0 = f(a) + f(b);
        x1 = 0;
        x2 = 0;
        for i = 1:n-1
            x = a + i*h;
            if bitget(i,1) x1 = x1 + f(x);
            else x2 = x2 + f(x);
            end
        end

        s = h*(x0 + 2*x2 + 4*x1)/3;
    end

    function s = trapezoid(f, a, b, n)
        h = (b-a)/n;
        x0 = f(a) + f(b);
        x1 = 0;
        for i = 1:n-1
            x = a + i*h;
            x1 = x1 + f(x);
        end

        s = h*(x0 + 2*x1)/2;
    end

    function s = midpoint(f, a, b, n)
        h = (b-a)/(n+2);
        x0 = 0;
        for i = 0:2:n
            x = a + (i+1)*h;
            x0 = x0 + f(x);
        end

        s = 2*h*x0;
    end

    s1 = trapezoid(f, a, b, n);
    s2 = simpson(f, a, b, n);
    s3 = midpoint(f, a, b, n);
end


f = @(x) x.^2.*log(x.^2+1);
a = 0; b = 2; n = 8;
[s1 s2 s3] = composite_integration(f, a, b, n);
disp(sprintf('Composite Trapezoidal rule: %1.9f', s1));
disp(sprintf('Composite Simpson rule: %1.9f', s2));
disp(sprintf('Composite Midpoint rule: %1.9f', s3));
