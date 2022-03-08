function [U, Nt, T, x, h, t] = diff_integration(L, T, Nx, tau, nu, u_0, f)
%% Secondary paramters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

Nt = ceil(T/tau) + 1;
T = tau*(Nt - 1);
t = 0 : tau : T;

%% Compact coefficients
alpha = 4*(6*nu + 5);
a = -2*(6*nu - 1)/alpha;
b = 4*(6*nu - 5)/alpha;
c = -2*(6*nu + 1)/alpha;

p = tau/alpha;
q = 10*tau/alpha;

%% Transition matrices
U_next = eye(Nx) + ...
    diag(a * ones(Nx-1, 1), 1) + ...
    diag(a * ones(Nx-1, 1),-1);
U_now = b*eye(Nx) + ...
    diag(c * ones(Nx-1, 1), 1) + ...
    diag(c * ones(Nx-1, 1),-1);

f_mat = q*eye(Nx) + ...
    diag(p * ones(Nx-1, 1), 1) + ...
    diag(p * ones(Nx-1, 1),-1);

%% Periodic border condition
U_next(1, end) = a;
U_next(end, 1) = a;

U_now(1, end) = c;
U_now(end, 1) = c;

f_mat(1, end) = p;
f_mat(end, 1) = p;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(x(1:end-1));

for k = 2 : Nt
    rp = -U_now*U(k - 1, :).' + ...
        f_mat*(f(t(k), x(1:end-1).') + f(t(k-1), x(1:end-1).'));
    U(k, :) = U_next \ rp;
end
U = [U, U(:, 1)];
end
