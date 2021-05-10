function [a b] = newton_divided_diff(x, f)
    %
end


x = [0.0; 0.2; 0.4; 0.6; 0.8];
f = [1.0000; 1.22140; 1.49182; 1.82212; 2.22554];
[a b] = newton_divided_diff(x, f);
disp(sprintf('a = [%1.5f %1.5f %1.5f %1.5f %1.5f]', a));
disp(sprintf('b = [%1.5f %1.5f %1.5f %1.5f %1.5f]', b));
