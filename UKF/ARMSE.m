function res = ARMSE(M,N,xj,xhj)
    mn = 1/(M*N);
    err = xj - xhj;
    res = sqrt(sum(err(:).^2 * mn));
endfunction
