function Xk = sigmap(x,P,c)
   A = c*Cholesky(P);
   Y = x(:, ones(1,numel(x)));
   Xk = [x Y+A Y-A];
endfunction
