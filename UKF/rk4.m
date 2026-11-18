function x_next = rk4(f_dot, x,u,dt)
    k1 = f_dot(x,u);
    k2 = f_dot(x + 0.5 * dt * k1,u);
    k3 = f_dot(x + 0.5 * dt * k2,u);
    k4 = f_dot(x + dt * k3,u);
    x_next = x + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
endfunction
