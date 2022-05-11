n1 = 400;
m1 = 4096;

D = 1;
t = 0.01/m1;
U = Schredinger(n1, m1, D, t);



os = zeros(1, 4);
n1 = 200;
m1 = 1024;
t = 0.01/m1;
for k = 1:4
    if (k >= 2)
        n1 = n1/2;
        m1 = m1/4;
    end

    t = 0.01/m1;
    U_test = Schredinger(n1, m1, D, t);
    
    Norm_eps = zeros(n1, m1+1);
    for i = 1:n1
        for j  =1:m1
            Norm_eps(i, j) = U_test(i, j)-U(2^k*i-(2^k-1), (4^k)*j-(4^k-1));
        end
   
    end

    maxx_eps = zeros(n1, 1);
    for i = 1:m1
        maxx_eps(i) = max(abs(Norm_eps(:, i)));
    end
    os(5-k) = max(maxx_eps);
    
end
Y = [25, 50, 100, 200];
 
loglog(Y, os);
title(['D =',num2str(D),', T =',num2str(t*m1), ',  nu = ', num2str(D*t*n1*n1)])
xlabel('N');
ylabel('log-error in C norm');
legend('Analytical');
