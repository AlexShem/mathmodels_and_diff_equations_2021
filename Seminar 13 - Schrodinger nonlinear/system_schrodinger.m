function [U, V] = system_schrodinger(u0, v0, D, h, tau, Nx, Nt, scheme)
U = zeros(Nt, Nx);
V = zeros(Nt, Nx);

U(1, :) = u0;
V(1, :) = v0;

for k = 2 : Nt
    if strcmp(scheme, 'euler')
        [u, v] = step_euler(U(k-1, :), V(k-1, :), D, h, tau);
        U(k, :) = u;
        V(k, :) = v;
    elseif strcmp(scheme, 'laksvendroff')
        error('not suppoted yet');
    elseif strcmp(scheme, 'adamsbashforth')
        [u, v] = step_adamsbashforth(U(k-1, :), V(k-1, :), D, h, tau);
        U(k, :) = u;
        V(k, :) = v;
    end
end
end
