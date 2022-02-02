clear
clc
syms a b c % Coeffs of the left side (approximating u)
syms p q r s % Coeffs of the right side (approximating f)
syms D h tau nu positive % Positive parameters
syms x t u(t, x) f(t, x)

coefs_u_app = [a b c];
coefs_f_app = [p q r s];

test_funs = [1, t, t^2, x^2, x^2*t, x^2*t^2, x^4];
test_funs = test_funs(:);



% Differential equation: Rod
du_fun = @(U) diff(U, t) - 1i*D*diff(U, x, 2);

coef_eqs = sym('Eqs', [length(test_funs), 1]);
% coef_eqs = zeros(11, 1);

for k = 1 : numel(test_funs)
    u(t, x) = test_funs(k);
    f(t, x) = du_fun(u(t, x));
    
    u_comact = ...
        u(t+tau, x) + a * (u(t+tau, x+h) + u(t+tau, x-h)) + ...
        b*u(t, x) + c * (u(t, x+h) + u(t, x-h));
    u_comact = subs(u_comact, [t, x], [0, 0]);
    
    f_comact = ...
        q*f(t+tau, x) + s * (f(t+tau, x+h) + f(t+tau, x-h)) + ...
        p*f(t, x) + r * (f(t, x+h) + f(t, x-h));
    f_comact = subs(f_comact, [t, x], [0, 0]);
    
    coef_eqs(k, 1) = u_comact == f_comact;
end

comp_eqs = solve(coef_eqs, [coefs_u_app, coefs_f_app]);

cfs = struct2array(comp_eqs).';
cfs = simplify(subs(cfs, D*tau, -1i*nu*h^2));

alpha = 24*nu + 20;
cfs_norm = simplify(cfs * alpha);
