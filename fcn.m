function [x]= fcn(u,v,uP,vP)
      x = []
      x(1,1) = uP;
      x(1,2) = vP;
      x(1,3) = -2*v/(1 + u*u + v*v)*(uP*vP) ;
      x(1,4) =vP + -2*u/(1 + u*u + v*v)*(uP*vP) ;
endfunction