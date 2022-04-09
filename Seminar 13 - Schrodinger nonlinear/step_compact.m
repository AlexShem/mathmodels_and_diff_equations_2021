function [U, V] = step_compact(u, v, up, vp, params)
D = params.D;
h = params.h;
tau = params.tau;

a = params.a;
b = 1;
c = (D*tau + 2*D*a*tau)/(2*h^2);
d = -(D*tau + 2*D*a*tau)/h^2;
p = -a*tau/2;
q = -tau/2;

C = params.C;
Dc = 1; % D coefficient (of the compact scheme)
A = -(D*tau + 2*D*C*tau)/(2*h^2);
B = -(D*tau + 2*D*C*tau)/h^2;
P = C*tau/2;
Q = tau/2;

F = -a*circshift(up, 1) - b*up - a*circshift(up, -1) ...
    -a*circshift(u, 1) - b*u - a*circshift(u, -1) ...
    -c*circshift(vp, 1) - d*vp - c*circshift(vp, -1) ...
    -c*circshift(vp, 1) - d*vp - c*circshift(vp, -1) + ...
    -p0*circshift(u, 1).^2/2 - q0*u.^2/2 - r0*circshift(u, -1).^2/2 ...
    -p1*circshift(u, 1).^2/2 - q1*u.^2/2 - r1*circshift(u, -1).^2/2;
A = diag(b1 + q1*u) + ...
    diag(c1 + r1*u(2:end), 1) + ...
    diag(a1 + p1*u(1:end-1), -1);
A(1, end) = a1 + p1*u(end);
A(end, 1) = c1 + r1*u(1);
Eps = A \ F.';
u = u + Eps.';

vp = v + tau/(2*h^2)*(circshift(u, 1) - 2*u + circshift(u, -1)) + ...
    D*tau/2 * u .* (u.^2 + v.^2);
up = u - tau/(2*h^2)*(circshift(v, 1) - 2*v + circshift(v, -1)) - ...
    D*tau/2 * v .* (u.^2 + v.^2);
vc = v + tau/h^2 * (circshift(up, 1) - 2*up + circshift(up, -1)) + ...
    D*tau * up .* (up.^2 + vp.^2);
uc = u - tau/h^2*(circshift(vp, 1) - 2*vp + circshift(vp, -1)) - ...
    D*tau * vp .* (up.^2 + vp.^2);
U = uc;
V = vc;
end
