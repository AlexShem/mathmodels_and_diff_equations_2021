%% Main parameters
C = 1;
L = 1;
T = 1;

Nx = 201;
nu = .1; % nu = C * tau / h

%% Secondary parameters
x = linspace(0, L, Nx);
h = x(2) - x(1);

tau = nu * h / C;
Nt = ceil(T/tau) + 1;

u_0 = sin(x*2*pi/L);
% sub_ind = x >= .45*L & x <= .55*L; 
% xsub = x(sub_ind);
% u_0 = zeros(size(x));
% u_0(sub_ind) = sin((xsub-L/2)*20/L*pi);
u_tau = u_0;

%% Transition matrices
U_next = eye(Nx);
U_now = (2-2*nu^2)*eye(Nx) + ...
    diag(nu^2 * ones(Nx-1, 1), 1) + ...
    diag(nu^2 * ones(Nx-1, 1),-1);
U_prev = -U_next;

%% Dirichlet border condition
% U_next(1, 2) = -1; U_next(end, end-1) = -1;
U_now(1, :) = 0;
U_now(end, :) = 0;
U_prev(1, :) = 0;
U_prev(end, :) = 0;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0;
U(2, :) = u_tau;

for k = 3 : Nt
    U(k, :) = U_next \ (U_now * U(k - 1, :).' + U_prev * U(k - 2, :).');
end

%% Visualisation
figure(2)

for k = [1 : floor(Nt/200) : Nt, Nt]
    plot(x, U(k, :));
    axis([0 L -1 1]);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end
