function [t,u,v,uP,vP] = difCen(a,b,y0,h)
  i=1;
  u(1) = y0(1,1);
  t(1)= a ;
  v(1) = y0(1,2);
  uP(1) = y0(1,3);
  vP(1) = y0(1,4);
  u(2) = u(1) + h*(uP(1));
  v(2) = v(1) + h*(vP(1));
  uP(2) = uP(1) + h*(-2*v(1)/(1 + u(1)*u(1) + v(1)*v(1)))*(uP(1)*vP(1)) ;
  vP(2) =vP(1) + h*(-2*u(1)/(1 + u(1)*u(1) + v(1)*v(1)))*(uP(1)*vP(1)) ;
  t(2) = t(1) + h;
  i = 2;
  f = ceil((b-a)/h);
  while i <= f+1
    i = i+1;
    t(i)= t(i-1)+h ;
    u(i) = u(i-2) + 2*h*(uP(i-1));
    v(i) = v(i-2) + 2*h*(vP(i-1));
    uP(i) = uP(i-2) + 2*h*(-2*v(i-1)/(1 + u(i-1)*u(i-1) + v(i-1)*v(i-1)))*(uP(i-1)*vP(i-1)) ;
    vP(i) =vP(i-2) + 2*h*(-2*u(i-1)/(1 + u(i-1)*u(i-1) + v(i-1)*v(i-1)))*(uP(i-1)*vP(i-1)) ;
  endwhile
endfunction
