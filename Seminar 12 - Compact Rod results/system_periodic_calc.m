function [U, T_fin] = system_periodic_calc(params, U_init, f, ...
    comp_scheme, full_sol, animate)
% Function that numericaly integrates the Rod Transverse Vibreation equation
%     params = [E, R, rho, T, L, Nx, tau]
%     U_init = structure with function handles
%         U   = @(x) value of U_0 at time t = 0
%         dU  = @(x) value of dU_0 at time t = 0
%     f == f(t, x) forcing; the right side
%     compute_scheme =  '555' for full compact scheme
%                       '353' for reduced compact sheme
%                       'CN'  for Crank-Nicolson scheme
%     full_sol = true if the full history is required
%              = false for final result only (improves memory)
%     animate  = true for the animation
%              = false for pure calculation

%% Extraction of the parameters
E = params(1);
R = params(2);
rho = params(3);
T = params(4);
L = params(5);
Nx = params(6);
tau = params(7);
U0 = U_init.U;
Utau = U_init.Utau;

h = L/Nx;
C = E*R^2/rho;
D = R^2;
nu = C*tau^2/h^4;
mu = D/h^2;

%% Chechs and definitions
switch comp_scheme
    case '555'
        if(nu <= 0 || mu < max(0, nu - 17/120))
            error('h and tau for compact scheme "555" do not satisfy  the stability constraint n < mu + 17/120');
        end
        alpha = 72*(41 + 30*nu + 90*mu);
        % Left part
        a = 96*(7 - 15*nu - 30*mu) / alpha;
        b = -144*(41 - 150*nu + 90*mu) / alpha;
        c = -192*(7 + 75*nu - 30*mu) / alpha;
        d = -24*(1 - 150*nu - 30*mu) / alpha;
        e = 12*(1 + 30*nu - 30*mu) / alpha;
        % Right part
        p0 = 10*tau^2 / alpha;
        q0 = 560*tau^2 / alpha;
        r0 = 2460*tau^2 / alpha;
        p1 = tau^2 / alpha;
        q1 = 56*tau^2 / alpha;
        r1 = 246*tau^2 / alpha;
    case '353'
        if(nu <= 0 || mu < max(0, nu - 1/12))
            error('h and tau for compact scheme "353" do not satisfy  the stability constraint n < mu + 10/12');
        end
        alpha = 96*(1 + 3*mu);
        % Left part
        a = 24*(1 - 6*mu) / alpha;
        b = -96*(2 - 9*nu + 6*mu) / alpha;
        c = -48*(1 + 12*nu - 6*mu) / alpha;
        d = 144*nu / alpha;
        e = 0;
        % Right part
        p0 = -10*tau^2*(nu - mu) / alpha;
        q0 = 20*tau^2*(1 + 2*nu - 2*mu) / alpha;
        r0 = 20*tau^2*(4 - 3*nu + 3*mu) / alpha;
        p1 = -tau^2*(nu - mu) / alpha;
        q1 = 2*tau^2*(1 + 2*nu - 2*mu) / alpha;
        r1 = 2*tau^2*(4 - 3*nu + 3*mu) / alpha;
    case 'CN'
        % (CN is absolutely stable)
        % Left part
        alpha = 1 + 3*nu + 2*mu;
        a = -(2*nu + mu) / alpha;
        b = -(2 + 4*mu) / alpha;
        c = 2*mu / alpha;
        d = 0;
        e = nu/2 / alpha;
        % Right part
        p0 = 0;
        p1 = 0;
        q0 = 0;
        q1 = 0;
        r0 = 0;
        r1 = tau^2/2 / alpha;
    otherwise
        error(['Undefined computing scheme ' comp_scheme]);
end

%% Initialisation of computing matrices;
U_next = sparse(diag(ones(Nx, 1)) + ...
    diag(a*ones(Nx - 1, 1), 1) + diag(a*ones(Nx - 1, 1), -1) + ...
    diag(e*ones(Nx - 2, 1), 2) + diag(e*ones(Nx - 2, 1), -2));
U_now = sparse(diag(b*ones(Nx, 1)) + ...
    diag(c*ones(Nx - 1, 1), 1) + diag(c*ones(Nx - 1, 1), -1) + ...
    diag(d*ones(Nx - 2, 1), 2) + diag(d*ones(Nx - 2, 1), -2));

f_next = sparse(diag(r1*ones(Nx, 1)) + ...
    diag(q1*ones(Nx - 1, 1), 1) + diag(q1*ones(Nx - 1, 1), -1) + ...
    diag(p1*ones(Nx - 2, 1), 2) + diag(p1*ones(Nx - 2, 1), -2));
f_now = sparse(diag(r0*ones(Nx, 1)) + ...
    diag(q0*ones(Nx - 1, 1), 1) + diag(q0*ones(Nx - 1, 1), -1) + ...
    diag(p0*ones(Nx - 2, 1), 2) + diag(p0*ones(Nx - 2, 1), -2));

% Left periodic condition
U_next(1, end-1:end) = [e a];
U_next(2, end) = e;
U_now(1, end-1:end) = [d c];
U_now(2, end) = d;
% Right periodic condition
U_next(end-1, 1) = e;
U_next(end, 1:2) = [a e];
U_now(end-1, 1) = d;
U_now(end, 1:2) = [c d];
U_prev = U_next;

% Left periodic condition
f_next(1, end-1:end) = [p1 q1];
f_next(2, end) = p1;
f_now(1, end-1:end) = [p0 q0];
f_now(2, end) = p0;
% Right periodic condition
f_next(end-1, 1) = p1;
f_next(end, 1:2) = [q1 p1];
f_now(end-1, 1) = p0;
f_now(end, 1:2) = [q0 p0];
f_prev = f_next;

%% Calculation of U
Nt = ceil(T/tau) + 1;
T_fin = tau*(Nt - 1);

x = linspace(0, L, Nx + 1); % +1: compensation for the L
x = x(1:end-1);             % (condiser every but the last point)
U0_val = U0(x);
% Utau_val = U0_val + tau * dU0(x);  % This can be improved
Utau_val = Utau(tau, x);
if full_sol
    U = zeros(Nt, Nx);
else
    U = zeros(3, Nx);
end

% Initial conditions
U(1, :) = U0_val;
U(2, :) = Utau_val;
f_val = zeros(3, Nx);
f_val(1, :) = f(0, x);
f_val(2, :) = f(tau, x);

if animate
    figure(10);
    plot([x, 2*x(end) - x(end-1)], [U(1, :), U(1, 1)]);
    axis([0, L, -1, 1]);
    title('t = 0 \cdot 10^{-3}');
end

for k = 3 : Nt
    f_val(3, :) = f(tau*(k-1), x);
    if full_sol
        right_part = -U_now*U(k-1, :).' - U_prev*U(k-2, :).' + ...
            f_next*f_val(3, :).' + f_now*f_val(2, :).' + f_prev*f_val(1, :).';
        U(k, :) = U_next \ right_part;
        
        if animate && ismember(k, [1:ceil(Nt/200):Nt, Nt])
            plot([x, 2*x(end) - x(end-1)], [U(k, :), U(k, 1)]);
            axis([0, L, -1, 1]);
            title(['t = ' num2str(tau*(k-1) * 1e3) '\cdot10^{-3}']);
            drawnow;
        end
    else
        right_part = -U_now*U(2, :).' - U_prev*U(1, :).' + ...
            f_next*f_val(3, :).' + f_now*f_val(2, :).' + f_prev*f_val(1, :).';
        U(3, :) = U_next \ right_part;
        U(1, :) = U(2, :);
        U(2, :) = U(3, :);
        
        if animate && ismember(k, [1:ceil(Nt/200):Nt, Nt])
            plot([x, 2*x(end) - x(end-1)], [U(3, :), U(3, 1)]);
            axis([0, L, -1, 1]);
            title(['t = ' num2str(tau*(k-1) * 1e3) '\cdot10^{-3}']);
            drawnow;
        end
    end
    f_val(1, :) = f_val(2, :);
    f_val(2, :) = f_val(3, :);
end

%% Returning the results
U = [U, U(:, 1)];
if ~full_sol
    U = U(end, :);
end

end
