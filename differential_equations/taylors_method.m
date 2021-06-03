function w = IVP_Taylor_Hermite(f, fp, a, b, alpha, N, t0)
    function T = Taylor(h, f_val, fp_val)
        T = f_val + h/2*fp_val;
    end
    
    function appr = Cubic_Hermite(t, t1, t2, h, w1, w2, f1, f2)
        c_0_1 = f1;
        c_1_2 = (w2 - w1) / h;
        c_2_3 = f2;
        
        c_0_1_2 = (c_1_2 - c_0_1) / h;
        c_1_2_3 = (c_2_3 - c_1_2) / h;
        
        c_0_1_2_3 = (c_1_2_3 - c_0_1_2) / h;
        appr = w1 + (t - t1)*c_0_1 + (t - t1)*(t - t1)*c_0_1_2 + (t - t1)*(t - t1)*(t - t2)*c_0_1_2_3;
    end
    
    h = (b-a) / N;
    w = alpha;
    t = a;
    for i=1:N
        f1 = f(t, w);
        fp1 = fp(t, w);
        T = Taylor(h, f1, fp1);
        if (t <= t0) && (t0 < t+h)
            w1 = w + h*T;
            t1 = t + h;
            fp2 = fp(t1, w1);
            f2 = f(t1, w1);
            break
        end
        w = w + h*T;
        t = t + h;
    end
    w = Cubic_Hermite(t0, t, t1, h, w, w1, f1, f2);
end


f = @(t,y) 2/t*y + t^2*exp(t);
fp = @(t,y) 2/t^2*y + (4*t + t^2)*exp(t);
a = 1;
b = 2;
alpha = 0;
N = 10;
t0 = 1.55;
w = IVP_Taylor_Hermite(f, fp, a, b, alpha, N, t0);
disp(sprintf('y(%1.4f): %2.7f', t0, w));
