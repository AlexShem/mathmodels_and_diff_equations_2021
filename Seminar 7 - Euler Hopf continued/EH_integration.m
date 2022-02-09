function [U, Nt, T, x, h, t] = EH_integration(L, T, Nx, tau, u_0, f, scheme)
%% Secondary paramters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

Nt = ceil(T/tau) + 1;
T = tau*(Nt - 1);
t = 0 : tau : T;

%% Compact coefficients
a1 = 1; b1 = 4; c1 = 1;
a0 = -1; b0 = -4; c0 = -1;

p1 = -3*tau/(2*h); q1 = 0; r1 = 3*tau/(2*h);
p0 = -3*tau/(2*h); q0 = 0; r0 = 3*tau/(2*h);

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(x(1:end-1));

for k = 2 : Nt
    u = U(k-1, :);
    U_LV = laks_vendroff(u, f, tau, h);
    if strcmp(scheme, 'compact')
        F = -a1*circshift(U_LV, 1) - b1*U_LV - c1*circshift(U_LV, -1) ...
            -a0*circshift(u, 1) - b0*u - c0*circshift(u, -1) ...
            -p0*circshift(u, 1).^2/2 - q0*u.^2/2 - r0*circshift(u, -1).^2/2 ...
            -p1*circshift(U_LV, 1).^2/2 - q1*U_LV.^2/2 - r1*circshift(U_LV, -1).^2/2;
        A = diag(b1 + q1*U_LV) + ...
            diag(c1 + r1*U_LV(2:end), 1) + ...
            diag(a1 + p1*U_LV(1:end-1), -1);
        A(1, end) = a1 + p1*U_LV(end);
        A(end, 1) = c1 + r1*U_LV(1);
        Eps = A \ F.';
        U_LV = U_LV + Eps.';
    end
    U(k, :) = U_LV;
end
U = [U, U(:, 1)];
end
