%% Function in ODE
% dx[theta(x) * dx u] + rho(x)*u = g(x)
theta = @(x) -ones(size(x));
dtheta = @(x) zeros(size(x));
rho = @(x) ones(size(x));
% If u(x) == exp(x), then g(x) is
u_true = @(x) exp(x);
g = @(x) zeros(size(x));

%% Parameters
a = 0;
b = 1;
L = b - a;
ua = u_true(a);
ub = u_true(b);

T = 1;
tau = .00001;
Nt = ceil(T/tau) + 1;

Nx = 101;
h = (b - a)/(Nx-1);
x = linspace(a, b, Nx);

nu = tau/h^2;
%% Transition matrices
U_now = (1 - tau - 2*nu) * eye(Nx) + ...
    diag(+nu * ones(Nx-1, 1), 1) + ...
    diag(+nu * ones(Nx-1, 1), -1);
U_now(1, :) = 0; U_now(end, :) = 0;
U_next = eye(Nx);

%% Integration
u_0 = (ub - ua)/(b-a) * x + ua;
u_0 = zeros(Nx, 1);
U = zeros(Nt, Nx);
U(1, :) = u_0;

for k = 2 : Nt
    right_part = U_now * U(k-1, :).';
    right_part(1) = ua;
    right_part(end) = ub;
    U(k, :) = U_next \ right_part;
end

%% Visualisation
figure(2)

for k = [1 : floor(Nt/200) : Nt, Nt]
    plot(x, U(k, :));
    hold on;
    plot(x, u_true(x), '--r');
    hold off;
    axis([0 L 0 3]);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end
