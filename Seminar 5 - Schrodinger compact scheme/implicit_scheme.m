%% Main parameters
D = 1;
L = 2*pi;
T = 1;

Nx = 10;
nu = .1i; % nu = D * tau / h^2

%% Secondary parameters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

tau = nu * h^2 / D / 1i;
Nt = ceil(T/tau) + 1;
T = tau*(Nt - 1);
t = 0 : tau : T;

u_ref = @(t, x) exp(-1i*D*t) .* cos(x);
u_0 = u_ref(0, x);
% figure(1)
% plot(x, u_0)

%% Transition matrices
alpha = 24*nu + 20;
a = 2*(1 - 6*nu)/alpha;
b = 4*(6*nu - 5)/alpha;
c = -2*(1 + 6*nu)/alpha;

U_now = b * eye(Nx) + ...
    diag(c * ones(Nx-1, 1), 1) + ...
    diag(c * ones(Nx-1, 1),-1);
U_next = eye(Nx) + ...
    diag(a * ones(Nx-1, 1), 1) + ...
    diag(a * ones(Nx-1, 1),-1);

%% Periodic border condition
U_next(1, end) = a;
U_next(end, 1) = a;

U_now(1, end) = c;
U_now(end, 1) = c;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(1 : end-1);

for k = 2 : Nt
    U(k, :) = -U_next \ U_now * U(k - 1, :).';
end
U = [U, U(:, 1)];

%% Visualisation
figure(2)
for k = [1 : floor(Nt/200) : Nt, Nt]
    plot(x, real(U(k, :)));
    hold on;
    plot(x, imag(U(k, :)));
    hold off;
    axis([0 L 0 1]);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end

%% Visualisation 2
prob = sqrt(h*sum(U .* conj(U), 2));

figure(3);
plot(t, prob);
xlabel('t');

%% Error
C_norm = max(abs(U - u_ref(t.', x)), [], 2);
figure(4)
plot(t, C_norm);
xlabel('t')
