function approx = newton_cotes(f, a, b, n, Closed)
    function approx = closed_newton_cotes(f, a, b, n)
        h = (b-a)/n;
        if n==1
            approx = h * (f(a) + f(b)) / 2;
        elseif n==2
            approx = h * (f(a) + 4*f(a+h) + f(b)) / 3;
        elseif n==3
            approx = 3*h * (f(a) + 3*f(a+h) + 3*f(a+2*h) + f(b)) / 8;
        else
            approx = 2*h * (7*f(a) + 32*f(a+h) + 12*f(a+2*h) + 32*f(a+3*h) + 7*f(b)) / 45;
        end
    end
    
    function approx = open_newton_cotes(f, a, b, n)
        h = (b-a)/(n+2);
        if n==0
            approx = 2*h * f(a+h);
        elseif n==1
            approx = 3*h * (f(a+h) + f(a+2*h)) / 2;
        elseif n==2
            approx = 4*h * (2*f(a+h) - f(a+2*h) + 2*f(a+3*h)) / 3;
        else
            approx = 5*h * (11*f(a+h) + f(a+2*h) + f(a+3*h) + 11*f(a+4*h)) / 24;
        end
    end

    if Closed
        approx = closed_newton_cotes(f, a, b, n);
    else
        approx = open_newton_cotes(f, a, b, n);
    end
end


f = @(x) (log(x))^3/3/x;
a = 2; b = 2.5; n = 3; Closed = 1;
approx = newton_cotes(f, a, b, n, Closed);
if Closed
    disp(sprintf('Closed Newton-Cotes (n=%d): %1.9f', n, approx));
else
    disp(sprintf('Open Newton-Cotes (n=%d): %1.9f', n, approx));
end
