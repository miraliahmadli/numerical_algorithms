format long;
% Setting up matrix A
[U,X]=qr(randn(80));
[V,X]=qr(randn(80));
S=exp(1).^(-1:-1:-80);
S=diag(S);
A=U*S*V;
 
% QR Factorization of A using 3 methods
[Q,RC]=cgs(A);
[Q,RM]=mgs(A);
[Q,RH]=qr(A);
 
% Visualization
figure(1);
loglog(diag(S),abs(diag(S)),'r-');
hold on;
loglog(diag(S),abs(diag(RC)),'b.');
loglog(diag(S),abs(diag(RM)),'r.');
loglog(diag(S),abs(diag(RH)),'k.');
legend('S','cgs','mgs', 'qs');
hold off

function [q, r] = cgs(a)
    [m,n] = size(a);
    q = zeros(m,n);
    r = zeros(n,n);
    for j = 1:n
        v = a(:,j);
        for i = 1:j-1
           r(i,j) = q(:, i)' * a(:, j);
           v = v - r(i,j)*q(:, i);
        end
        r(j,j) = norm(v);
        if r(j,j) ~= 0
            q(:, j) = v / r(j,j);
        end
    end
end

function [q, r] = mgs(a)
    [m,n] = size(a);
    v = zeros(m,n);
    q = zeros(m,n);
    r = zeros(n,n);
    for i = 1:n
        v(:, i) = a(:, i);
    end
    for i = 1:n
        r(i,i) = norm(v(:, i));
        if r(i,i) ~= 0
            q(:, i) = v(:, i) / r(i,i);
        end
        for j = i+1:n
           r(i,j) = q(:,i)' * v(:,j);
           v(:,j) = v(:,j) - r(i,j)*q(:,i);
        end
    end
end
