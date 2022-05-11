firstnode = @(x, y)[y(2); 5*y(1)-5*y(1)^3]; % функция в компактном виде

xspan = [0, 2.3];
y0 = [1/2, 0];
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[x, y] = ode45(firstnode, xspan, y0, opts);

%plot(y(:, 1), y(:, 2))

y1 = zeros(74, 2);
x1 = zeros(74, 1);
for k = 1:74
    y1(k, 1) = y(10*k, 1);
    y1(k, 2) = y(10*k, 2);
    x1(k) = x(10*k);
end

plot(y1(:, 1), y1(:, 2))

xlabel('v');
ylabel('dv/dx')
