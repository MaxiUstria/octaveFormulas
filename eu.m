function [t,u,v,uP,vP] = eu(a,b,y0,h)
  i=1;
  u(1) = y0(1,1);
  t(1)= a ;
  v(1) = y0(1,2);
  uP(1) = y0(1,3);
  vP(1) = y0(1,4);
  while t(i)<b+h
    i=i+1;
    t(i)= t(i-1)+h ;
    u(i) = u(i-1) + h*(uP(i-1));
    v(i) = v(i-1) + h*(vP(i-1));
    uP(i) = uP(i-1) + h*(-2*v(i-1)/(1 + u(i-1)*u(i-1) + v(i-1)*v(i-1)))*(uP(i-1)*vP(i-1)) ;
    vP(i) = vP(i-1) + h*(-2*u(i-1)/(1 + u(i-1)*u(i-1) + v(i-1)*v(i-1)))*(uP(i-1)*vP(i-1)) ;
  endwhile
endfunction
