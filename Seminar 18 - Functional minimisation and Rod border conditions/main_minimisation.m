%% Main parameters
A = 0; % Left border value
B = 0; % Right border value

theta = @(x) ones(size(x));
f = @(x) sin(x);

F = @(x, u, du) theta(x).*du.^2/2 + f(x).*u(x);
