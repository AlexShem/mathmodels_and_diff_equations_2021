function [U, Nt, T, x, h, t] = B_integration(L, T, D, Nx, tau, nu, u_0, scheme)
%% Secondary paramters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

Nt = ceil(T/tau) + 1;
T = tau*(Nt - 1);
t = 0 : tau : T;

%% Compact coefficients
a1 = -2/3 + 2*nu; b1 = -8/3 - 4*nu; c1 = -2/3 + 2*nu;
a0 = 2/3 + 2*nu; b0 = 8/3 - 4*nu; c0 = 2/3 + 2*nu;

p1 = tau/h; q1 = 0; r1 = -tau/h;
p0 = tau/h; q0 = 0; r0 = -tau/h;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(x(1 : end-1));

for k = 2 : Nt
    u = U(k-1, :);
    U_AB = adams_bashforth(u, D, tau, h);
    if strcmp(scheme, 'compact')
        F = -a1*circshift(U_AB, 1) - b1*U_AB - c1*circshift(U_AB, -1) ...
            -a0*circshift(u, 1) - b0*u - c0*circshift(u, -1) ...
            -p0*circshift(u, 1).^2/2 - q0*u.^2/2 - r0*circshift(u, -1).^2/2 ...
            -p1*circshift(U_AB, 1).^2/2 - q1*U_AB.^2/2 - r1*circshift(U_AB, -1).^2/2;
        A = diag(b1 + q1*U_AB) + ...
            diag(c1 + r1*U_AB(2:end), 1) + ...
            diag(a1 + p1*U_AB(1:end-1), -1);
        A(1, end) = a1 + p1*U_AB(end);
        A(end, 1) = c1 + r1*U_AB(1);
        Eps = A \ F.';
        U_AB = U_AB + Eps.';
    end
    U(k, :) = U_AB;
end
U = [U, U(:, 1)];
end
