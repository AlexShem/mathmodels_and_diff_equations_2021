%% Parameters

scheme = 'compact'; params.r = 0.1;
% scheme = 'standard';

L = 2*pi;

Nx = 100;
Ny = Nx;

hx = L / Nx;
hy = L / Ny;

x = linspace(0, L, Nx);
y = linspace(0, L, Ny);
[x, y] = meshgrid(x, y);

params.x = x;
params.y = y;
params.hx = hx;
params.hy = hy; 
params.Nx = Nx;
params.Ny = Ny;

%% Right-hand side
syms xs ys
Lu = @(u) diff(u, xs, 2) + diff(u, ys, 2);

% u_sym = sin(xs)*sin(ys);
u_sym = exp(-(xs-L/2)^2 - (ys-L/2)^2);

f_sym = Lu(u_sym);
u_true = matlabFunction(u_sym, 'Vars', [xs, ys]);
f = matlabFunction(f_sym, 'Vars', [xs, ys]);
clear xs ys u_sym f_sym

%% Integration
u = system_poisson_dirichlet(scheme, params, f);

%% Visualisation
figure(1)
% contour(x, y, u, 'ShowText', 'on');
surf(x, y, u);
xl = xlabel('$x$'); xl.Interpreter = 'latex'; xl.FontSize = 16;
yl = ylabel('$y$'); yl.Interpreter = 'latex'; yl.FontSize = 16;

%% Visualisation 2
figure(2)
contour(x, y, log10(abs(u - u_true(x, y))), 18, 'ShowText', 'on');
xl = xlabel('$x$'); xl.Interpreter = 'latex'; xl.FontSize = 16;
yl = ylabel('$y$'); yl.Interpreter = 'latex'; yl.FontSize = 16;
ttl = title(['$\log_{10} |u - u^*|$']); ttl.Interpreter = 'latex'; ttl.FontSize = 16;
