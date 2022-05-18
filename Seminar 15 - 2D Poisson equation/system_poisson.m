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
    a = -1/5;
    b = -1/20;
    c = 1;
    r = params.r;
    p = -hx^2/5 + 4*r;
    q = -hx^2/40 - 2*r;
else
    error(['Scheme "', scheme, '" is not supported yet.']);
end

A = c*eye(N);
F = p*eye(N);

%% X Left border
% X Left border: first
A(1, 1 + Ny) = a; % right
A(1, 2) = a; % down
A(1, 1 + N - Ny) = a; %left
A(1, Ny) = a; %up
A(1, 1 + 2*Ny - 1) = b; % tr
A(1, 1 + Ny + 1) = b; % br
A(1, N - Ny + 2) = b; % bl
A(1, N) = b; % tl

F(1, 1 + Ny) = q; % right
F(1, 2) = q; % down
F(1, 1 + N - Ny) = q; %left
F(1, Ny) = q; %up

% X Left border: intermediate
for k = 2 : Ny-1
    % a coeffs
    A(k, k + Ny) = a; % right
    A(k, k + 1) = a; % down
    A(k, k + N - Ny) = a; %left
    A(k, k - 1) = a; %up
    A(k, k + Ny - 1) = b; % tr
    A(k, k + Ny + 1) = b; % br
    A(k, k + N - Ny + 1) = b; % bl
    A(k, k + N - Ny - 1) = b; % tl
    
    F(k, k + Ny) = q; % right
    F(k, k + 1) = q; % down
    F(k, k + N - Ny) = q; %left
    F(k, k - 1) = q; %up
end

% X Left border: last
A(Ny, 2*Ny) = a; % right
A(Ny, 1) = a; % down
A(Ny, Ny + N - Ny) = a; %left
A(Ny, Ny - 1) = a; %up
A(Ny, Ny + Ny - 1) = b; % tr
A(Ny, Ny + 1) = b; % br
A(Ny, N - Ny + 1) = b; % bl
A(Ny, N - 1) = b; % tl

F(Ny, 2*Ny) = q; % right
F(Ny, 1) = q; % down
F(Ny, Ny + N - Ny) = q; %left
F(Ny, Ny - 1) = q; %up

%% Intermediate X steps
for k = Ny+1 : N-Ny
    if mod(k, Ny) == 1 % upper border
        A(k, k + Ny) = a; % right
        A(k, k + 1) = a; % down
        A(k, k - Ny) = a; %left
        A(k, k + Ny - 1) = a; %up
        A(k, k + 2*Ny - 1) = b; % tr
        A(k, k + Ny + 1) = b; % br
        A(k, k - Ny + 1) = b; % bl
        A(k, k - 1) = b; % tl
        
        F(k, k + Ny) = q; % right
        F(k, k + 1) = q; % down
        F(k, k - Ny) = q; %left
        F(k, k + Ny - 1) = q; %up
    elseif mod(k, Ny) == 0 % bottom border
        A(k, k + Ny) = a; % right
        A(k, k - Ny + 1) = a; % down
        A(k, k - Ny) = a; %left
        A(k, k - 1) = a; %up
        A(k, k + Ny - 1) = b; % tr
        A(k, k + 1) = b; % br
        A(k, k - 2*Ny + 1) = b; % bl
        A(k, k - Ny - 1) = b; % tl
        
        F(k, k + Ny) = q; % right
        F(k, k - Ny + 1) = q; % down
        F(k, k - Ny) = q; %left
        F(k, k - 1) = q; %up
    else % interior point
        A(k, k + Ny) = a; % right
        A(k, k + 1) = a; % down
        A(k, k - Ny) = a; %left
        A(k, k - 1) = a; %up
        A(k, k + Ny - 1) = b; % tr
        A(k, k + Ny + 1) = b; % br
        A(k, k - Ny + 1) = b; % bl
        A(k, k - Ny - 1) = b; % tl
        
        F(k, k + Ny) = q; % right
        F(k, k + 1) = q; % down
        F(k, k - Ny) = q; %left
        F(k, k - 1) = q; %up
    end
end

%% X Right border
% X Right border: first
A(N - Ny + 1, 1) = a; % right
A(N - Ny + 1, N - Ny + 1 + 1) = a; % down
A(N - Ny + 1, N - Ny + 1 - Ny) = a; %left
A(N - Ny + 1, N) = a; %up
A(N - Ny + 1, Ny) = b; % tr
A(N - Ny + 1, 2) = b; % br
A(N - Ny + 1, N - Ny + 1 - Ny + 1) = b; % bl
A(N - Ny + 1, N - Ny) = b; % tl

F(N - Ny + 1, 1) = q; % right
F(N - Ny + 1, N - Ny + 1 + 1) = q; % down
F(N - Ny + 1, N - Ny + 1 - Ny) = q; %left
F(N - Ny + 1, N) = q; %up

% X Right border: intermediate
for k = N - Ny+2 : N-1
    % a coeffs
    A(k, mod(k, Ny)) = a; % right
    A(k, k + 1) = a; % down
    A(k, k - Ny) = a; %left
    A(k, k - 1) = a; %up
    A(k, mod(k, Ny) - 1) = b; % tr
    A(k, mod(k, Ny) + 1) = b; % br
    A(k, k - Ny + 1) = b; % bl
    A(k, k - Ny - 1) = b; % tl
    
    F(k, mod(k, Ny)) = q; % right
    F(k, k + 1) = q; % down
    F(k, k - Ny) = q; %left
    F(k, k - 1) = q; %up
end

% X Right border: last
A(N, Ny) = a; % right
A(N, N - Ny + 1) = a; % down
A(N, N - Ny) = a; %left
A(N, N - 1) = a; %up
A(N, Ny - 1) = b; % tr
A(N, 1) = b; % br
A(N, N - 2*Ny + 1) = b; % bl
A(N, N - Ny - 1) = b; % tl

F(N, Ny) = q; % right
F(N, N - Ny + 1) = q; % down
F(N, N - Ny) = q; %left
F(N, N - 1) = q; %up

%% Solution
rhs = f(params.x, params.y);
u = A \ rhs(:);
u = reshape(u, Ny, Nx);
end
