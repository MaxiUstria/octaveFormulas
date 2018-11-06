function octv= fv(t,x)
      octv(1,1) = x(3);
      octv(1,2) = x(4);
      octv(1,3) = ((-2*x(2)/(1 + x(1)*x(1) + x(2)*x(2)))*(x(3)*x(4))) ;
      octv(1,4) = ((-2*x(1)/(1 + x(1)*x(1) + x(2)*x(2)))*(x(3)*x(4)));
      octv = octv(:)';
endfunction