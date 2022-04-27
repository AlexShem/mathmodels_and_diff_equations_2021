function [U, V] = step_maccormack(u, v, params)
D = params.D;
h = params.h;
tau = params.tau;

vp = v + tau/(h^2)*(circshift(u, 1) - 2*u + circshift(u, -1)) + ...
    tau * u .* (u.^2 + v.^2);
up = u - tau/(h^2)*(circshift(v, 1) - 2*v + circshift(v, -1)) - ...
    tau * v .* (u.^2 + v.^2);
vc = (v + vp)/2 + tau/(2*h^2) * (circshift(up, 1) - 2*up + circshift(up, -1)) + ...
    tau/2 * up .* (up.^2 + vp.^2);
uc = (u + up)/2 - tau/(2*h^2)*(circshift(vp, 1) - 2*vp + circshift(vp, -1)) - ...
    tau/2 * vp .* (up.^2 + vp.^2);
U = uc;
V = vc;
end
