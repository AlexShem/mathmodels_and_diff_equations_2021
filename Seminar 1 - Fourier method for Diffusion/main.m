L = 2;
D = .2;
K = 10;
T = 1;
a = 1; b = 2;

x = linspace(0, L, 101);
t = linspace(0, T, 201);
[x, t] = meshgrid(x, t);

u_0 = @(x) sin(2*pi/L * x);
C = zeros(K, 1);
e_k = @(x, k) sin(2*pi*k/L*x);

u = zeros(size(x));

for k = 1 : K
    fe = integral(@(x) u_0(x) .* e_k(x, k), 0, L);
    ee = integral(@(x) e_k(x, k).^2, 0, L);
    C(k) = fe / ee;
    u = u + C(k) * exp(-D * (2*pi*k/L)^2 * t) .* sin(2*pi*k/L * x);
    u_true = u + a + (b - a)/L * x;
end

figure(1);
surf(x, t, u, 'EdgeColor', 'none');
xlabel('x');
ylabel('t');

figure(2);
for j = 1 : size(t, 1)
    plot(x(j, :), u_true(j, :));
    title(['t = ', num2str(t(j, 1))]);
    axis([0, L, -1, 3])
    drawnow limitrate nocallbacks;
end
