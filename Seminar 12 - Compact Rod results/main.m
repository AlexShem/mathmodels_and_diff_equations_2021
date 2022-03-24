%% Set of problems
T = 1;

% Define the paramaters of the problems
L = 2*pi;
E = 210e9;
R = 1e-2;
rho = 7900;

C = E*R^2/rho;   % C = E*R^2/rho
D = R^2;

% Function to find
V = sqrt(C/(D + 1));
U_true_fn = @(t, x) sin(x + V*t);
U0_true_fn = @(x) sin(x);
f_true_fn = @(t, x) 0*t;
% U0_true_fn = @(x) t.^2.*sin(x);
% f_true_fn = @(t, x) (E + C*t.^2).*sin(x);

U_init = struct('U', @(x) U_true_fn(0, x), ...
    'Utau', @(tau, x) U_true_fn(tau, x));

full_sol = true;
animate = false;

%% Error in cycle
compute_scheme = {'CN', 'Compact_555'};

% C_norm_dt = cell(1, 2);
% L2_norm_dt = cell(1, 2);

mu_range = linspace(.01, .15, 5);
nu_range = linspace(.01, .15, 6);
[nu_mesh, mu_mesh] = meshgrid(nu_range, mu_range);

h_mesh = sqrt(D./mu_mesh);
Nx_mesh = ceil(L ./ h_mesh);
h_mesh = L ./ Nx_mesh;
tau_mesh = sqrt(nu_mesh .* h_mesh.^4 / C);

mu_mesh = D./h_mesh.^2;
nu_mesh = C * tau_mesh.^2 ./ h_mesh.^4;

C_norm_dt = struct(...
    'CN', zeros(length(mu_range), length(nu_range)),...
    'CN_time', zeros(length(mu_range), length(nu_range)),...
    'Compact_555', zeros(length(mu_range), length(nu_range)),...
    'Compact_555_time', zeros(length(mu_range), length(nu_range))...
    );

L2_norm_dt = struct(...
    'CN', zeros(length(mu_range), length(nu_range)),...
    'CN_time', zeros(length(mu_range), length(nu_range)),...
    'Compact_555', zeros(length(mu_range), length(nu_range)),...
    'Compact_555_time', zeros(length(mu_range), length(nu_range))...
    );

for comp_ind = 1:numel(compute_scheme)
    comp_sch = compute_scheme{comp_ind};
    tStart = tic;
    
    for i = 1 : length(mu_range) % Interations over rows
        for j = 1 : length(nu_range)
            Nx = Nx_mesh(i, j);
            tau = tau_mesh(i, j);
            
            x = linspace(0, L, Nx + 1);
            h = h_mesh(i, j);
            
            mu = D/h^2;
            nu = C*tau^2/h^4;
            
            params = [E, R, rho, T, L, Nx, tau];
            
            %%% Integration
            try
                [U, T_fin] = system_periodic_calc(params, U_init, f_true_fn, ...
                    replace(comp_sch, 'Compact_', ''), full_sol, animate);
                %%% Error study
                t = 0 : tau : T_fin;
                Nt = length(t);
                [x_mesh, t_mesh] = meshgrid(x, t);
                
                [Cn, Ct_ind] = max(max(abs(U - U_true_fn(t_mesh, x_mesh)), [], 2));
                C_norm_dt.(comp_sch)(i, j) = Cn;
                C_norm_dt.([comp_sch, '_time'])(i, j) = t(Ct_ind);
                
                [L2, L2t_ind] = max(sqrt(h/L * sum((U - U_true_fn(t_mesh, x_mesh)).^2, 2)));
                L2_norm_dt.(comp_sch)(i, j) = L2;
                L2_norm_dt.([comp_sch, '_time'])(i, j) = t(L2t_ind);
            catch
                C_norm_dt.(comp_sch)(i, j) = NaN;
                C_norm_dt.([comp_sch, '_time'])(i, j) = NaN;
                L2_norm_dt.(comp_sch)(i, j) = NaN;
                L2_norm_dt.([comp_sch, '_time'])(i, j) = NaN;
            end
        end
    end
    tEnd = toc(tStart);
    disp(['[Done]: ', comp_sch, ...
        '. Time ellapsed: ', num2str(tEnd)]);
end

% Save the results
filename = ['C_L2_' num2str(length(nu_range)) 'x' num2str(length(mu_range)) '.mat'];
C_norm = C_norm_dt;
L2_norm = L2_norm_dt;
save(filename, 'C_norm', 'L2_norm', ...
    'nu_mesh', 'mu_mesh', 'tau_mesh', 'Nx_mesh', 'h_mesh');

%% Load the results
load C_L2_6x5.mat;

%% Visualize the results
figure(1)
subplot(1, 2, 1);
contour(nu_mesh, mu_mesh, log10(C_norm.CN), 'ShowText', 'on');
xlabel('\nu'); ylabel('\mu');
title('CN')

subplot(1, 2, 2);
contour(nu_mesh, mu_mesh, log10(C_norm.Compact_555), 'ShowText', 'on');
xlabel('\nu'); ylabel('\mu');
title('5-5-5')
