%% Main parameters
f = @(x, y) x.*exp(-x.^2-y.^2)+(x.^2+y.^2)/20; % Function to minimise
fx = @(x) x(1).*exp(-x(1).^2-x(2).^2)+(x(1).^2+x(2).^2)/20;
lb = [-2, -2];
ub = [2, 2];

figure(1)
fsurf(f, [lb(1) ub(1) lb(2) ub(2)], 'ShowContours', 'on', 'FaceAlpha', .3)
xlabel('x');
ylabel('y')

%% Gradient descent
x0 = [1.5; 1.5];
% x0 = [0; 0];
[xmin_g, fmin_g, niter_g, path_g] = grad_descent(x0, fx, [], 1000);

figure(1)
hold on;
plot3(path_g(1,:), path_g(2,:), f(path_g(1,:), path_g(2,:)), '-*r')
hold off;

%% Matlab native
option = optimoptions('fminunc', 'Display', 'iter', 'PlotFcn', 'optimplotfval');
[x, fval] = fminunc(fx, x0, option);

