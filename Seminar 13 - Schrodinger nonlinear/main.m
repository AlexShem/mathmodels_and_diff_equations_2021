%% Parameters

% scheme = 'compact';
% scheme = 'laksvendroff';
scheme = 'adamsbashforth';
% scheme = 'euler';

L = 2*pi;
T = 2*pi;

% nu = tau/h^2;
nu = 0.06;

Nx = 100;
% Nt = 2^12;
D = 0.01;

h = L / Nx;
% tau = T / Nt;
tau = nu*h^2;
Nt = ceil(T / tau) + 1;
T = tau*(Nt - 1);

x = linspace(0, L, Nx + 1);
x = x(1 : end-1);
t = linspace(0, T, Nt);

%% Initial conditions
u0 = sin(x);
v0 = zeros(size(x));

%% Integration
[u, v] = system_schrodinger(u0, v0, D, h, tau, Nx, Nt, scheme);
I = sum(u.^2 + v.^2, 2);

%% Visualisation
figure(1)
plot(t, I)
xlabel('t');
