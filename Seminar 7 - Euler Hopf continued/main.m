%% Main parameters
L = 2*pi;
T = .5;
nu = 0.05; % nu = tau/h
Nx = [25, 50, 100, 200];
scheme = 'compact';

u_0 = @(x) sin(x);
f = @(u) u.^2/2;

Nx_ref = 400;
h_ref = L / (Nx_ref - 1);
tau_ref = nu * h_ref;
[U_ref, ~, ~, x_ref, ~, t_ref] = EH_integration(L, T, Nx_ref, tau_ref, u_0, f, scheme);
[x_ref, t_ref] = meshgrid(x_ref, t_ref);

%%
C_norm = zeros(size(Nx));
figure(1);

for j = 1 : length(Nx)
    x = linspace(0, L, Nx(j) + 1);
    h = x(2) - x(1);
    tau = nu * h;
    
%     [U, ~, ~, x, ~, t] = EH_integration(L, T, Nx(j), tau, u_0, f, scheme);
    [U, ~, ~, x, ~, t] = EH_smooth(L, T, Nx(j), tau, u_0, f, scheme);
    [x_mesh, t_mesh] = meshgrid(x, t);
    U_ref_val = interp2(x_ref, t_ref, U_ref, x_mesh, t_mesh);
    Cn = max(abs(U_ref_val - U), [], 2);
    plot(t, Cn); hold on;
    
    Cn = Cn(~isnan(Cn));
    C_norm(j) = Cn(end);
end
hold off;

%%
figure(2);
loglog(Nx, C_norm);

%%
order_comp = (log10(C_norm(end)) - log10(C_norm(1))) / ...
    (log10(Nx(end)) - log10(Nx(1)));
disp(['[ORDER] ', scheme, ' scheme: ', num2str(order_comp), '. nu = ', num2str(nu)])
