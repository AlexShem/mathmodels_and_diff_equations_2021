function u = system_poisson(scheme, params, f)
hx = params.hx;
% hy = params.hy;
N = numel(params.x);
Nx = params.Nx;
Ny = params.Ny;

if strcmp(scheme, 'standard')
    a = 1;
    b = 0;
    c = -4;
    p = hx^2;
    q = 0;
    r = 0;
elseif strcmp(scheme, 'compact')
else
    error(['Scheme "', scheme, '" is not supported yet.']);
end

A = c*eye(N);
F = p*eye(N);

%% X Left border
% X Left border: first
A(1, 1 + Ny) = a; % right
A(1, 2) = a; % down
A(1, 1 + (Nx-1)*Ny) = a; %left
A(1, Ny) = a; %up

% X Left border: intermediate
for k = 2 : Ny-1
    % a coeffs
    A(k, k + Ny) = a; % right
    A(k, k + 1) = a; % down
    A(k, k + (Nx-1)*Ny) = a; %left
    A(k, k - 1) = a; %up
end

% X Left border: last
A(Ny, 2*Ny) = a; % right
A(Ny, 1) = a; % down
A(Ny, Ny + (Nx-1)*Ny) = a; %left
A(Ny, Ny - 1) = a; %up

%% Intermediate X steps
for k = Ny+1 : (Nx-1)*Ny
    if mod(k, Ny) == 1 % upper border
        A(k, k + Ny) = a; % right
        A(k, k + 1) = a; % down
        A(k, k - Ny) = a; %left
        A(k, k + Ny - 1) = a; %up
    elseif mod(k, Ny) == 0 % bottom border
        A(k, k + Ny) = a; % right
        A(k, k - Ny + 1) = a; % down
        A(k, k - Ny) = a; %left
        A(k, k - 1) = a; %up
    else % interior point
        A(k, k + Ny) = a; % right
        A(k, k + 1) = a; % down
        A(k, k - Ny) = a; %left
        A(k, k - 1) = a; %up
    end
end

%% X Right border
% X Right border: first
A((Nx-1)*Ny + 1, 1) = a; % right
A((Nx-1)*Ny + 1, (Nx-1)*Ny + 1 + 1) = a; % down
A((Nx-1)*Ny + 1, (Nx-1)*Ny + 1 - Ny) = a; %left
A((Nx-1)*Ny + 1, N) = a; %up

% X Right border: intermediate
for k = (Nx-1)*Ny+2 : N-1
    % a coeffs
    A(k, mod(k, Ny)) = a; % right
    A(k, k + 1) = a; % down
    A(k, k - Ny) = a; %left
    A(k, k - 1) = a; %up
end

% X Right border: last
A(N, Ny) = a; % right
A(N, N - Ny + 1) = a; % down
A(N, N - Ny) = a; %left
A(N, N - 1) = a; %up
end
