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

N=400;
M=60;

f_dot=@(x,u)[
          alp*x(2)*x(3) + u(1)/J1;
          bet*x(1)*x(3) + u(2)/J2;
          gam*x(2)*x(1) + u(3)/J3
      ];

h=@(x)[x(1); x(2);x(3)]


xVr = zeros(n,N,M);
sVr = zeros(n,N,M);

for j = 1:M
  s=[0.3;0.8;1];
  x=s+q*randn(n,1);
  P = eye(n);
  xV = zeros(n,N);
  sV = zeros(n,N);
  zV = zeros(m,N);
  for k=1:N
    tc = (k-1)*dt;

    u = [0.2*sin(2*tc); -0.1*cos(tc); 0.15];

    f = @(x)rk4(f_dot,x,u,dt);
    z = h(s) + r*randn(m,1);
    sV(:,k)= s;
    zV(:,k) = z;

    [x, P] = ukf(f,x,P,h,z,Q,R);
    xV(:,k) = x;
    s = f(s) + q*randn(n,1);
  endfor
  xVr(:,:,j) = xV;
  sVr(:,:,j) = sV;
endfor

arm = ARMSE(M,N,sVr,xVr)

for k=1:3
subplot(3,1,k)
plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
legend('Истинные значения','UKF')
end
