function x_next = rk4(f_dot, x, dt)
    k1 = f_dot(x);
    k2 = f_dot(x + 0.5 * dt * k1);
    k3 = f_dot(x + 0.5 * dt * k2);
    k4 = f_dot(x + dt * k3);
    x_next = x + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
endfunction
