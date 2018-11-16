tic
h = [0.1 0.05 0.025 0.01 0.005 0.0025 0.001 0.0005 0.00025 0.0001 ];
y0 = [-9 -8 3 2];
t0 = 0;
tf = 2;
for i = 1 : 10 
  N = (tf-t0)/h(i);
  [t,backward] = beul("oct","dfcn",y0,t0,tf,N);
  [t,trapecio] = trap("oct","dfcn",y0,t0,tf,N);
  [t,heunn] = heun("oct",t0,tf,y0,N);
  [t,rk44] = rk4("oct",t0,tf,y0,N);
  [t,solver] = ode45("fv", t, y0);
  [t u v uP vP] = eu(t0,tf,y0,h(i));
  [t,rk22] = rk2("oct",t0,tf,y0,N);
  eulerr = [];
  eulerr(1:N+1,1) = u(1:N+1);
  eulerr(1:N+1,2) = v(1:N+1);
  eulerr(1:N+1,3) = uP(1:N+1);
  eulerr(1:N+1,4) = vP(1:N+1);
  errore(i) = mean(mean(abs(eulerr - solver)));
  errorb(i) = mean(mean(abs(backward - solver)));
  errort(i) = mean(mean(abs(trapecio - solver)));
  errorh(i) = mean(mean(abs(heunn - solver)));
  errorrk(i) = mean(mean(abs(rk44 - solver)));
  errordc(i) = mean(mean(abs(rk22 - solver)));
  errore2(i) = (norm(eulerr - solver,inf));
  errorb2(i) = (norm(backward - solver,inf));
  errort2(i) = (norm(trapecio - solver, inf));
  errorh2(i) = (norm(heunn - solver,inf));
  errorrk2(i) = (norm(rk44 - solver,inf));
  errordc2(i) = (norm(rk22 - solver,inf));
end
toc