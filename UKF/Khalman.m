clear;
n=3;
m=3;

q=0.01;
r=0.1;
Q=q^2*eye(n);
R=r^2*eye(m);

dt = 0.01;

J1 = 0.5;
J2 = 0.3;
J3 = 0.1;

alp = -(J3-J2)/J1;
bet = -(J1-J3)/J2;
gam = -(J2-J1)/J3;

N=200;
M=50;

f_dot=@(x,u)[
          alp*x(2)*x(3) + u(1)/J1;
          bet*x(1)*x(3) + u(2)/J2;
          gam*x(2)*x(1) + u(3)/J3
      ];

h=@(x)[x(1); x(2);x(3)]

w = [sqrt(3)/4; -sqrt(3); 1/sqrt(3)];

w1e = w(1);
w2e = w(2);
w3e = w(3);

ue = -[J1 * alp * w2e * w3e;
       J2 * bet * w1e * w3e;
       J3 * gam * w1e * w2e];

A = [0,       alp*w3e,  alp*w2e;
     bet*w1e, 0,        bet*w3e;
     gam*w2e, gam*w1e,  0
    ];

B = diag([1/J1,1/J2,1/J3]);

Hk = eye(n);

EM = expm([A,B; zeros(m,n+m)]*dt);
Fk = EM(1:n,1:n);
Bk = EM(1:n,n+1:end);

xVr = zeros(n,N,M);
sVr = zeros(n,N,M);

for j = 1:M
  s=[0.3;0.8;1];
  x=s+q*randn(n,1);
  xk = x - w;
  Pk = 0.1*eye(n);
  xV = zeros(n,N);
  sV = zeros(n,N);
  zV = zeros(m,N);
  for k=1:N
    u = [0.2; 0.1; 0.15];
    mu = u - ue;
    f = @(x)rk4(f_dot,x,u,dt);
    z = h(s) + r*randn(m,1);
    sV(:,k)= s;
    zV(:,k) = z;

    xk = Fk * xk + Bk * mu;
    Pk = Fk * Pk * Fk.' + Q;

    yk = (z-w) - Hk*xk;
    Sk = Hk*Pk*Hk.' + R;
    Kk = Pk*Hk.' / Sk;

    xk = xk + Kk*yk;
    Pk = (eye(n) - Kk*Hk)*Pk;

    x = w + xk;
    xV(:,k) = x;

    s = f(s) +q*randn(n,1);

  endfor
  xVr(:,:,j) = xV;
  sVr(:,:,j) = sV;
endfor

arm = ARMSE(M,N,sVr,xVr)
titles = {"Движение спутника; Уравнения Эйлера \n w(1) - угловая скорость", "w(2)", "w(3)"};
for k=1:3
subplot(3,1,k)
plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--', 1:N, zV(k,:), '---')
title(titles(k));
legend('Истинные значения','UKF','Зашумлённые данные')
end



