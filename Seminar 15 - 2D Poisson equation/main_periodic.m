%% Parameters

% scheme = 'compact'; params.a = -.51; params.C = -.51;
scheme = 'standard';

L = 2*pi;
% T = 12;

Nx = 4;
Ny = Nx;

hx = L / Nx;
hy = L / Ny;

x = linspace(0, L, Nx + 1);
x = x(1 : end-1);
y = linspace(0, L, Ny + 1);
y = y(1 : end-1);
[x, y] = meshgrid(x, y);

params.x = x;
params.y = y;
params.hx = hx;
params.hy = hy;
params.Nx = Nx;
params.Ny = Ny;

%% Right-hand side
f = @(x, y) zeros(size(x));
% f = @(x, y) -2*sin(x).*sin(y);

%% Integration
u = system_poisson(scheme, params, f);
x = [x, L];
u = [u, u(:, 1)];
v = [v, v(:, 1)];
I = sum(u.^2 + v.^2, 2);

%% Visualisation
figure(1)
hold on;
plot(t, I)
xl = xlabel('$t$'); xl.Interpreter = 'latex'; xl.FontSize = 16;
yl = ylabel('$I(u, v) = u^2 + v^2$'); yl.Interpreter = 'latex'; yl.FontSize = 16;
hold off;
ttl = title(['$\tau = $', num2str(tau), ', $h = $', num2str(hx)]); ttl.Interpreter = 'latex'; ttl.FontSize = 16;

%% Visualisation 2
figure(2)
hold on;
plot(t, log10(I./I(1)))
xl = xlabel('$t$'); xl.Interpreter = 'latex'; xl.FontSize = 16;
yl = ylabel('$\log_{10} I(u, v) / I_0(u, v)$'); yl.Interpreter = 'latex'; yl.FontSize = 16;
hold off;
ttl = title(['$\tau = $', num2str(tau), ', $h = $', num2str(hx)]); ttl.Interpreter = 'latex'; ttl.FontSize = 16;

%% Animation
figure(2)
for k = [1 : 100 : Nt, Nt]
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
