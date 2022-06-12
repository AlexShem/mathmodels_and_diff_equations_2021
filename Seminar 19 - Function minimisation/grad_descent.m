function [fmin, niter, path] = grad_descent(x0, f, tol, maxiter)
if nargin < 3
    tol = 1e-6;
end
if nargin < 4
    maxiter = 50;
end

x0 = x0(:);
path = zeros(length(x0), maxiter + 1);
path(:, 1) = x0;

grad_x = grad_f(x0, f);
x = x0 - grad_x;
iter = 1;
path(:, iter + 1) = x;

df = abs(f(x) - f(x0));

while (df > tol) && iter < maxiter
    grad_x = grad_f(x, f);
    grad_x0 = grad_f(x0, f);
    lam = abs((x - x0).' * (grad_x - grad_x0)) / norm(grad_x - grad_x0, 2);

    x0 = x;
    x = x0 - lam*grad_x;
    iter = iter + 1;
    df = abs(f(x) - f(x0));
    path(:, iter + 1) = x;
end

if iter == maxiter && df > tol
    warning('Gradient descent algorithm did not converge properly');
end
fmin = f(x);
niter = iter;
path = path(:, 1 : niter + 1);

function grad = grad_f(x, f)
d = length(x);
h = 1e-8;
grad = arrayfun(@(ind) (f(x + h*((1:d).' == ind)) - f(x))/h, (1:d).');
