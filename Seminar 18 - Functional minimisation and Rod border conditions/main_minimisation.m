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

%% Minimisation problem
% rng('default');
u0 = -1 + 2*rand(length(x), 1);
u0 = zeros(length(x), 1);
u0(1) = uA; u0(end) = uB;
Aeq = zeros(2, Nx); Aeq(1, 1) = 1; Aeq(2, Nx) = 1;
beq = [uA; uB];

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set', 'PlotFcn', 'optimplotfval', 'MaxFunctionEvaluations', 1000*Nx, 'TolFun', 1e-8);
% options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
% options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set');

tic
[u_min, Fval, exitflag, output] = fmincon(@(u) objectivefun(x, F, u), u0, [], [], Aeq, beq, [], [], [], options);
toc

%% Functional
function Fval = objectivefun(x, F, u)
h = x(2) - x(1);
du = (u(3:end) - u(1:end-2)) / (2*h);
du = [(u(2)-u(1))/h; du; (u(end)-u(end-1))/h];

Fval = F(x, u, du);
Fval = h*sum(Fval);
end
