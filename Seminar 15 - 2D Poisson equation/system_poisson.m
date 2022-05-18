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
A(Ny, 1 + Ny) = a; % right
A(Ny, 1) = a; % down
A(Ny, 1 + (Nx-1)*Ny) = a; %left
A(Ny, Ny) = a; %up
end
