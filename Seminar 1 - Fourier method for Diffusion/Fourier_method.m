D = .2;
L = 1;
T = 1;
K = 10; % Количество членов в ряде Фурье
A = .5; B = 1;

x = linspace(0, L, 11);
t = linspace(0, T, 101);
[x, t] = meshgrid(x, t);

u_0 = @(x) sin(pi/L * x);
% u_0 = @(x) x.*(L - x);
C = zeros(K, 1);

u = zeros(size(x));
for k = 1 : K
    fe = integral(@(x) u_0(x) .* sin(pi * k * x / L), 0, L);
    ee = integral(@(x) sin(pi * k * x / L).^2, 0, L);
    C(k) = fe / ee;
    
    u = u + C(k) * exp(-D * (pi*k/L)^2 * t) .* sin(pi*k*x/L);
end

u_AB = u + A + (B - A)/L * x;

figure(1);
surf(x, t, u_AB, 'EdgeColor', 'none', 'FaceColor', 'interp');
xlabel('x');
ylabel('t');

figure(2);
for j = 1 : size(t, 1)
    plot(x(j, :), u_AB(j, :));
    axis([0 L 0 2]);
    xlabel('x');
    title(['t = ', num2str(t(j, 1))]);
    drawnow;
%     pause(.1)
end
