%% Parameters

scheme = 'compact'; params.a = -1; params.C = -1;
% scheme = 'laksvendroff';
% scheme = 'maccormack';
% scheme = 'adamsbashforth';
% scheme = 'euler';

L = 2*pi;
T = .2*pi;

Nx = 5;
nu = .004; % nu = tau/h^2
D = 1;

h = L / Nx;
tau = nu*h^2;
Nt = ceil(T / tau) + 1;
T = tau*(Nt - 1);

x = linspace(0, L, Nx + 1);
x = x(1 : end-1);
t = linspace(0, T, Nt);

params.D = D;
params.h = h;
params.tau = tau;
params.Nx = Nx;
params.Nt = Nt;

%% Initial conditions
u0 = sin(x);
% v0 = cos(x);
v0 = zeros(size(x));

%% Integration
[u, v] = system_schrodinger(u0, v0, params, scheme);
x = [x, L];
u = [u, u(:, 1)];
v = [v, v(:, 1)];
I = sum(u.^2 + v.^2, 2);

%% Visualisation
figure(1)
hold on;
plot(t, I)
xlabel('t');
hold off;

%% Animation
figure(2)
for k = [1 : 10 : Nt, Nt]
    plot(x, u(k, :));
    hold on;
    plot(x, v(k, :));
    hold off;
    lg = legend('$u$', '$v$', 'Interpreter', 'latex');
    lg.FontSize = 16; lg.Location = 'southeast';
    title(['Scheme: ', scheme]);
    subtitle(['t = ', num2str((k-1)*tau)]);
    xlabel('$x$, m.', 'Interpreter', 'latex', 'FontSize', 14);
    axis([min(x), max(x), -1.1, 1.1])
    
    drawnow;
end
