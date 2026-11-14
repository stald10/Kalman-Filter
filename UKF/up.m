function [x,Xk,P,X1] = up(f,X,Wm,Wc,n,R)
  L =  size(X,2);
  x =  zeros(n,1);
  Xk = zeros(n,L);
  for k=1:L
    Xk(:,k)=f(X(:,k));
    x = x + Wm(k)*Xk(:,k);
  endfor
  X1 = Xk - x(:,ones(1,L));
  P = X1*diag(Wc)*X1.' + R
endfunction
