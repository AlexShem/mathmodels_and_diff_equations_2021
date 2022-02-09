function P = desmooth(u)
    P = u - (circshift(u, -1) + circshift(u, 1))/2;
end
