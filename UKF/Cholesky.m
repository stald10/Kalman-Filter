function L = Cholesky(P)
      n = size(P,1);
      T = zeros(n,n);
      for i = 1:n
            s = 0;
            for p = 1:i-1
              s = s + T(i,p)^2;
            endfor
            T(i,i) = sqrt(P(i,i) - s);
        for j = i+1:n
            s = 0;
            for p = 1:i-1
              s = s + T(i,p)*T(j,p);
            endfor
            T(j,i) = (1/T(i,i))* (P(j,i) - s);
        endfor
      endfor
      L = T;
endfunction
