function [a b] = divided_diff_for_Hermite(x, f, g)
    %
end


x = [8.0; 8.6];
f = [17.56492; 18.50515];
g = [3.116256; 3.151762];
[a b] = divided_diff_for_Hermite(x, f, g);
disp(sprintf('a = [%1.5f %1.5f %1.5f %1.5f]', a));
disp(sprintf('b = [%1.5f %1.5f %1.5f %1.5f]', b));
