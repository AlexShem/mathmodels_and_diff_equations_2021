%% Main parameters
D = .5i;
L = 2*pi;
T = 1;
test_type = 'homogeneous';

Nx = 25:25:200;
nu = .5i; % nu = D * tau / h^2

switch test_type
    case 'homogeneous'
        % Homogeneous case
        u_ref = @(t, x) exp(-D*t) .* cos(x);
        u_0 = @(x) u_ref(0, x);
        f = @(t, x) zeros(max(size(x), 1));
    case 'nonhomogeneous'
        % Nonhomogeneous case
        u_ref = @(t, x) cos(x) .* sin(t);
        u_0 = @(x) u_ref(0, x);
        f = @(t, x) cos(x).*(cos(t) + D*sin(t));
    otherwise
        error(['Problem "' test_type '" is not supported.']);
end

%% Integration
C_norm = zeros(size(Nx));
figure(1);
if ~isreal(D) && strcmp(test_type, 'homogeneous')
    figure(2);
end

for j = 1 : length(Nx)
    x = linspace(0, L, Nx(j) + 1);
    h = x(2) - x(1);
    tau = nu * h^2 / D;
    
    [U, ~, ~, x, ~, t] = diff_integration(L, T, Nx(j), tau, nu, u_0, f);
    [x_mesh, t_mesh] = meshgrid(x, t);
    U_ref_val = u_ref(t_mesh, x_mesh);
    Cn = max(abs(U_ref_val - U), [], 2);
    
    figure(1)
    semilogy(t, Cn); hold on;
    
    Cn = Cn(~isnan(Cn));
    C_norm(j) = Cn(end);
    
    if ~isreal(D) && strcmp(test_type, 'homogeneous')
        prob = sqrt(h*sum(abs(U).^2, 2));
        figure(2);
        plot(t, prob); hold on;
    end
end

figure(1); hold off; legend(num2str(Nx.'))
if ~isreal(D) && strcmp(test_type, 'homogeneous')
    figure(2); hold off; legend(num2str(Nx.'))
end

%%
figure(3);
loglog(Nx, C_norm, '-o');
xlabel('N_x');
ylabel('C norm of |u_{ref} - u|')

%%
order_comp = (log10(C_norm(end)) - log10(C_norm(1))) / ...
    (log10(Nx(end)) - log10(Nx(1)));
disp(['[ORDER] ', num2str(-order_comp), '. nu = ', num2str(nu)])
