Fk = [1.0,1.0; 0.0, 1.0];
Bk = [0.0; 1.0];

Qk = 0.001*eye(2);

Rk = 25.0;

Hk = [1.0,0.0];


mu = 0;
sigma = 5;

truth = 100 + 3.0*(0:59);

nk = normrnd(mu,sigma,1,length(truth));

zk = truth + nk;

xk = [zk(1);0.0];

uk = 0.0;

Pk = [25.0,0; 0.0, 100.0];

err = [xk(1)];

for measurement  = zk(2:end)
  xk = Fk * xk + Bk * uk;
  Pk = Fk * Pk * Fk.' + Qk;

  yk = measurement  - Hk*xk;

  Sk = Hk*Pk*Hk.' + Rk;

  Kk = Pk*Hk.'* inv(Sk);

  xk = xk + Kk*yk;

  Pk = (eye(2) - Kk*Hk)*Pk;

  err = [err,xk(1)];

endfor
fprintf("Detection error: %.2f px\n", sqrt(mean((zk - truth) .^ 2)));
fprintf("Kalman error: %.2f px\n", sqrt(mean((err - truth) .^ 2)));
fprintf("Estimated velocity: %.2f px/frame\n", xk(2));


tk = 1/6*(0:59);

plot(tk,truth,"r-;Истинные данные;",tk, zk, "g-;Зашумлённые данные;",tk,err, "b-;Офильтрованные данные;");



