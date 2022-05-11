function [U, V] = step_compact(u, v, up, vp, params)
h = params.h;
tau = params.tau;
Nx = params.Nx;

a = params.a;
b = 1;
c = (tau + 2*a*tau)/(2*h^2);
e = -(tau + 2*a*tau)/h^2;
p = a*tau/2;
q = tau/2;

C = params.C;
E = 1; % E coefficient (of the compact scheme)
A = -(tau + 2*C*tau)/(2*h^2);
B = (tau + 2*C*tau)/h^2;
P = C*tau/2;
Q = tau/2;

u0 = circshift(u, 1);
u2 = circshift(u, -1);
v0 = circshift(v, 1);
v2 = circshift(v, -1);
up0 = circshift(up, 1);
up2 = circshift(up, -1);
vp0 = circshift(vp, 1);
vp2 = circshift(vp, -1);

F = -a*(up0 + up2) - b*up ...
    +a*(u0 + u2) + b*u ...
    -c*(vp0 + vp2) - e*vp ...
    -c*(v0 + v2) - e*v ...
    -p*(up0.^2.*vp0 + vp0.^3 + up2.^2.*vp2 + vp2.^3) - q*(up.^2.*vp + vp.^3) ...
    -p*(u0.^2.*v0 + v0.^3 + u2.^2.*v2 + v2.^3) - q*(u.^2.*v + v.^3);
G = -A*(up0 + up2) - B*up ...
    -A*(u0 + u2) - B*u ...
    -C*(vp0 + vp2) - E*vp ...
    +C*(v0 + v2) + E*v ...
    +P*(vp0.^2.*up0 + up0.^3 + vp2.^2.*up2 + up2.^3) + Q*(vp.^2.*up + up.^3) ...
    +P*(v0.^2.*u0 + u0.^3 + v2.^2.*u2 + u2.^3) + Q*(v.^2.*u + u.^3);
rhs = [F, G].';

% Epsilon Delta matrices
eps_11 = diag(b + 2*q*up.*vp) ...
    + diag(a + 2*p*up(2:end).*vp(2:end), 1) ...
    + diag(a + 2*p*up(1:end-1).*vp(1:end-1), -1);
delta_12 = diag(e + q*(up.^2 + 3*vp.^2)) ...
    + diag(c + p*(up(2:end).^2 + 3*vp(2:end).^2), 1) ...
    + diag(c + p*(up(1:end-1).^2 + 3*vp(1:end-1).^2), -1);
eps_21 = diag(B - Q*(vp.^2 + 3*up.^2)) ...
    + diag(A - P*(vp(2:end).^2 + 3*up(2:end).^2), 1) ...
    + diag(A - P*(vp(1:end-1).^2 + 3*up(1:end-1).^2), -1);
delta_22 = diag(E - 2*Q*up.*vp) ...
    + diag(C - 2*P*up(2:end).*vp(2:end), 1) ...
    + diag(C - 2*P*up(1:end-1).*vp(1:end-1), -1);

% Border conditions for the matrices
eps_11(1, end) = a + 2*p*up(end).*vp(end);
eps_11(end, 1) = a + 2*p*up(1).*vp(1);
delta_12(1, end) = c + p*(up(end).^2 + 3*vp(end).^2);
delta_12(end, 1) = c + p*(up(1).^2 + 3*vp(1).^2);
eps_21(1, end) = A - P*(vp(end).^2 + 3*up(end).^2);
eps_21(end, 1) = A - P*(vp(1).^2 + 3*up(1).^2);
delta_22(1, end) = C - 2*P*up(end).*vp(end);
delta_22(end, 1) = C - 2*P*up(1).*vp(1);

ED_mat = [eps_11, delta_12;
    eps_21, delta_22];
ED = ED_mat \ rhs;

epsilon = ED(1:Nx);
delta = ED(Nx+1:end);

U = up + epsilon.';
V = vp + delta.';
end
