function [U, V] = step_maccormack(u, v, params)
h = params.h;
tau = params.tau;

% Calculate the predicted values
up = u - tau/h^2*(circshift(v, 1) - 2*v + circshift(v, -1)) - tau*v.*(u.^2 + v.^2);
vp = v + tau/h^2*(circshift(u, 1) - 2*u + circshift(u, -1)) + tau*u.*(u.^2 + v.^2);

U = .5*(u + up) - .5*tau/h^2*(circshift(vp, 1) - 2*vp + circshift(vp, -1)) - .5*tau*vp.*(up.^2 + vp.^2);
V = .5*(v + vp) + .5*tau/h^2*(circshift(up, 1) - 2*up + circshift(up, -1)) + .5*tau*up.*(up.^2 + vp.^2);
end
