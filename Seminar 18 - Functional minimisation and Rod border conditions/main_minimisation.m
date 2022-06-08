%% Main parameters
uA = 0; % Left border value
uB = 0; % Right border value
L = 2*pi;

theta = @(x) ones(size(x));
% theta = @(x) x;
f = @(x) sin(x);
u_c = @(x, c1, c2) c1 + c2.*x - sin(x);

F = @(x, u, du) theta(x).*du.^2/2 + f(x).*u;

Nx = 105;
x = linspace(0, L, Nx);
x = x(:);
h = x(2) - x(1);

%% Divergent scheme
xh = .5*(x(2:end) + x(1:end-1));
th = theta(xh);

D = -diag([0; th] + [th; 0]) + ...
    diag(th, 1) + ...
    diag(th, -1);
D(1, :) = 0; D(1, 1) = 1;
D(end, :) = 0; D(end, end) = 1;

rhs = h^2*f(x);
rhs(1) = uA; rhs(end) = uB;

u_ds = D \ rhs;
du_ds = (u_ds(3:end) - u_ds(1:end-2)) / (2*h);
du_ds = [(u_ds(2)-u_ds(1))/h; du_ds; (u_ds(end)-u_ds(end-1))/h];
F_ds = h*sum(F(x, u_ds, du_ds));
