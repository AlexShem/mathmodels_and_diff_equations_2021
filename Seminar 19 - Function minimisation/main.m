%% Main parameters

f = @(x, y) x.*exp(-x.^2-y.^2)+(x.^2+y.^2)/20; % Function to minimise
figure(1)
fsurf(f, [-2, 2], 'ShowContours', 'on', 'FaceAlpha', .3)
xlabel('x');
ylabel('y')

fx = @(x) x(1).*exp(-x(1).^2-x(2).^2)+(x(1).^2+x(2).^2)/20;

%% Gradient descent
x0 = [1.5; 1.5];
[xmin, fmin, niter, path] = grad_descent(x0, fx, [], 1000)

figure(1)
hold on;
plot3(path(1,:), path(2,:), f(path(1,:), path(2,:)), '-*r')
hold off;

