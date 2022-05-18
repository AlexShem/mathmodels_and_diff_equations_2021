function u = system_poisson_dirichlet(scheme, params, f)
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
A(1, 1 + Ny + 1) = b; % br

F(1, 1 + Ny) = q; % right
F(1, 2) = q; % down
F(1, 1 + Ny + 1) = r; % br

% X Left border: intermediate
for k = 2 : Ny-1
    A(k, k + Ny) = a; % right
    A(k, k + 1) = a; % down
    A(k, k - 1) = a; %up
    A(k, k + Ny - 1) = b; % tr
    A(k, k + Ny + 1) = b; % br
    
    F(k, k + Ny) = q; % right
    F(k, k + 1) = q; % down
    F(k, k - 1) = q; %up
    F(k, k + Ny - 1) = r; % tr
    F(k, k + Ny + 1) = r; % br
end

% X Left border: last
A(Ny, 2*Ny) = a; % right
A(Ny, Ny - 1) = a; %up
A(Ny, Ny + Ny - 1) = b; % tr

F(Ny, 2*Ny) = q; % right
F(Ny, Ny - 1) = q; %up
F(Ny, Ny + Ny - 1) = r; % tr

%% Intermediate X steps
for k = Ny+1 : N-Ny
    if mod(k, Ny) == 1 % upper border
        A(k, k + Ny) = a; % right
        A(k, k + 1) = a; % down
        A(k, k - Ny) = a; %left
        A(k, k + Ny + 1) = b; % br
        A(k, k - Ny + 1) = b; % bl
        
        F(k, k + Ny) = q; % right
        F(k, k + 1) = q; % down
        F(k, k - Ny) = q; %left
        F(k, k + Ny + 1) = r; % br
        F(k, k - Ny + 1) = r; % bl
    elseif mod(k, Ny) == 0 % bottom border
        A(k, k + Ny) = a; % right
        A(k, k - Ny) = a; %left
        A(k, k - 1) = a; %up
        A(k, k + Ny - 1) = b; % tr
        A(k, k - Ny - 1) = b; % tl
        
        F(k, k + Ny) = q; % right
        F(k, k - Ny) = q; %left
        F(k, k - 1) = q; %up
        F(k, k + Ny - 1) = r; % tr
        F(k, k - Ny - 1) = r; % tl
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
        F(k, k + Ny - 1) = r; % tr
        F(k, k + Ny + 1) = r; % br
        F(k, k - Ny + 1) = r; % bl
        F(k, k - Ny - 1) = r; % tl
    end
end

%% X Right border
% X Right border: first
A(N - Ny + 1, N - Ny + 1 + 1) = a; % down
A(N - Ny + 1, N - Ny + 1 - Ny) = a; %left
A(N - Ny + 1, N - Ny + 1 - Ny + 1) = b; % bl

F(N - Ny + 1, N - Ny + 1 + 1) = q; % down
F(N - Ny + 1, N - Ny + 1 - Ny) = q; %left
F(N - Ny + 1, N - Ny + 1 - Ny + 1) = r; % bl

% X Right border: intermediate
for k = N - Ny+2 : N-1
    A(k, k + 1) = a; % down
    A(k, k - Ny) = a; %left
    A(k, k - 1) = a; %up
    A(k, k - Ny + 1) = b; % bl
    A(k, k - Ny - 1) = b; % tl
    
    F(k, k + 1) = q; % down
    F(k, k - Ny) = q; %left
    F(k, k - 1) = q; %up
    F(k, k - Ny + 1) = r; % bl
    F(k, k - Ny - 1) = r; % tl
end

% X Right border: last
A(N, N - Ny) = a; %left
A(N, N - 1) = a; %up
A(N, N - Ny - 1) = b; % tl

F(N, N - Ny) = q; %left
F(N, N - 1) = q; %up
F(N, N - Ny - 1) = r; % tl

%% Solution
f_val = f(params.x, params.y);
rhs = F*f_val(:);
u = A \ rhs;
u = reshape(u, Ny, Nx);
end
