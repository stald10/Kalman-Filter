function [x,P]=ukf(f,x,P,h,z,Q,R)
  L = numel(x);
  m = numel(z);
  alpha = 1e-3;
  ki = 0;
  beta = 2;
  lambda = alpha^2*(L+ki)-L;
  c = L + lambda;
  Wm = [lambda/c 0.5/c + zeros(1,2*L)];
  Wc = Wm;
  Wc(1) = Wc(1) + (1-alpha^2 + beta);
  c = sqrt(c);
  X = sigmap(x,P,c);

  [x1,X1,P1,X2] = up(f,X,Wm,Wc,L,Q);
  [z1,Z1,P2,Z2] = up(h,X1,Wm,Wc,m,R);
  P12 = X2*diag(Wc)*Z2.';
  K = P12*inv(P2);
  x = x1 + K*(z-z1);
  P = P1 - K*P12';
endfunction

