%% Main parameters
uA = 0; % Left border value
uB = 0; % Right border value
L = 2*pi;

theta = @(x) ones(size(x));
f = @(x) sin(x);

F = @(x, u, du) theta(x).*du.^2/2 + f(x).*u(x);

Nx = 105;
x = linspace(0, L, Nx);
h = x(2) - x(1);
