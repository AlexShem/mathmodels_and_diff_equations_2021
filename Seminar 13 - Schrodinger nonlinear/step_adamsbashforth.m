function [U, V] = step_adamsbashforth(u, v, params)
D = params.D;
h = params.h;
tau = params.tau;

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
