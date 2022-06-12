%% Main parameters

f = @(x, y) x.*exp(-x.^2-y.^2)+(x.^2+y.^2)/20; % Function to minimise
figure(1)
fsurf(f, [-2, 2], 'ShowContours', 'on')
xlabel('x');
ylabel('y')

fx = @(x) x(1).*exp(-x(1).^2-x(2).^2)+(x(1).^2+x(2).^2)/20;

%% Gradient descent
x0 = [1; 1];
[fmin, niter, path] = grad_descent(x0, fx)

figure(1)
hold on;
plot3(path(1,:), path(2,:), f(path(1,:), path(2,:)), '-*r')
hold off;
