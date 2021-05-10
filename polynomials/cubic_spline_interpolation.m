function [a b c d] = cubic_spline(x, f, opt, g0, gn)
    %
end


x = [0; 1; 2; 3];
f = [1; exp(1); exp(2); exp(3)];
disp('===Natural cubic spline===');
[a b c d] = cubic_spline(x, f, 0);
disp(sprintf('a = [%1.5f %1.5f %1.5f]', a));
disp(sprintf('b = [%1.5f %1.5f %1.5f]', b));
disp(sprintf('c = [%1.5f %1.5f %1.5f]', c));
disp(sprintf('d = [%1.5f %1.5f %1.5f]', d));
