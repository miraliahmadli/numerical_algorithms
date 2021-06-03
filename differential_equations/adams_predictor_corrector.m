function w = IVP_Adams_Predictor_Corrector(f, a, b, alpha, N)
    h = (b-a) / N;
    w0 = alpha;
    t = a;
    w = zeros(N,1);
    w(1)=w0;
    
    % Runge-Kutta
    ws = zeros(1,4);
    ws(1) = w0;
    ts = zeros(1,4);
    ts(1) = a;
    for i=1:min([3, N])
        K1 = h*f(ts(i),ws(i)); 
        K2 = h*f(ts(i) + h/2, ws(i) + K1/2);
        K3 = h*f(ts(i) + h/2, ws(i) + K2/2);
        K4 = h*f(ts(i) + h, ws(i) + K3);
        ws(i+1) = ws(i) + (K1 + 2*K2 + 2*K3 + K4)/6;
        ts(i+1) = t+i*h;
        w(i+1) = ws(i+1);
    end
    
    for i=4:N
        t = a + i*h;
        v = 55*f(ts(4), ws(4)) - 59*f(ts(3), ws(3)) +37*f(ts(2), ws(2)) - 9*f(ts(1), ws(1));
        w0 = ws(4) + h*v/24;
        v = 9*f(t, w0) + 19*f(ts(4), ws(4)) - 5*f(ts(3), ws(3)) + f(ts(2), ws(2));
        w0 = ws(4) + h*v/24;
        
        for j=1:3
            ts(j) = ts(j+1);
            ws(j) = ws(j+1);
        end
        ts(4) = t;
        ws(4) = w0;
        w(i+1) = w0;
    end
end


f = @(t,y) (2-2*t*y)/(t^2+1);
a = 0;
b = 1;
alpha = 1;
N = 10;
w = IVP_Adams_Predictor_Corrector(f, a, b, alpha, N);
disp(sprintf('y(%1.4f): %2.7f', b, w(end)));
