function [U, V] = system_schrodinger(u0, v0, params, scheme)
Nx = params.Nx;
Nt = params.Nt;

U = zeros(Nt, Nx);
V = zeros(Nt, Nx);

U(1, :) = u0;
V(1, :) = v0;

for k = 2 : Nt
    if strcmp(scheme, 'euler')
        [u, v] = step_euler(U(k-1, :), V(k-1, :), params);
        U(k, :) = u;
        V(k, :) = v;
    elseif strcmp(scheme, 'laksvendroff')
        error(['Scheme ' scheme ' is not suppoted yet']);
    elseif strcmp(scheme, 'adamsbashforth')
        [u, v] = step_adamsbashforth(U(k-1, :), V(k-1, :), params);
        U(k, :) = u;
        V(k, :) = v;
    elseif strcmp(scheme, 'maccormack')
        [u, v] = step_maccormack(U(k-1, :), V(k-1, :), params);
        U(k, :) = u;
        V(k, :) = v;
    elseif strcmp(scheme, 'compact')
%         [up, vp] = step_adamsbashforth(U(k-1, :), V(k-1, :), params);
        [up, vp] = step_maccormack(U(k-1, :), V(k-1, :), params);
        [u, v] = step_compact(U(k-1, :), V(k-1, :), up, vp, params);
        U(k, :) = u;
        V(k, :) = v;
    end
end
end
