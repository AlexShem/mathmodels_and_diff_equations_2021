function [U, T, tau] = wave_integration(T, C, L, Nx, nu, u_0, u_tau)
% Secondary parameters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);
tau = nu * h / C;
Nt = ceil(T/tau) + 1;
T = tau*(Nt - 1);

% Compact coefficients
a = -(nu^2-1) / (2*(nu^2+5));
b = 10*(nu^2-1) / (nu^2+5);
c = -(5*nu^2+1) / (nu^2+5);

% Transition matrices
U_next = eye(Nx) + ...
    diag(a * ones(Nx-1, 1), 1) + ...
    diag(a * ones(Nx-1, 1),-1);
U_now = b * eye(Nx) + ...
    diag(c * ones(Nx-1, 1), 1) + ...
    diag(c * ones(Nx-1, 1),-1);

% Periodic border condition
U_next(1, end) = a;
U_next(end, 1) = a;
U_now(1, end) = c;
U_now(end, 1) = c;

U_prev = U_next;

% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(x(1 : end-1));
U(2, :) = u_tau(x(1 : end-1), tau);

for k = 3 : Nt
    U(k, :) = -U_next \ (U_now * U(k - 1, :).' + U_prev * U(k - 2, :).');
end
U = [U, U(:, 1)];
end

