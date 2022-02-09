%% Main parameters
L = 2*pi;
T = .99;
nu = 0.1; % nu = tau/h
Nx = 200;
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);
tau = nu * h;

u_0 = @(x) sin(x);
f = @(u) u.^2/2;

ax = [0 L -2 2];

%%
[U_ref, ~, ~, x_ref, ~, t_ref] = EH_integration(L, T, 400, tau, u_0, f, 'compact');

%% Visualisation
figure(2)
for k = [1 : floor(Nt/200) : Nt, Nt]
    plot(x, U(k, :));
    axis(ax);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end
