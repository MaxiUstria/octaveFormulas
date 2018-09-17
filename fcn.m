function oct= f(x,t);
      oct(1) = x(1);
      oct(2) = x(2);
      oct(3) = -2*x(2)/(1 + x(1)*x(1) + x(2)*x(2))*(x(3)*x(4)) ;
      oct(4) = -2*x(1)/(1 + x(1)*x(1) + x(2)*x(2))*(x(3)*x(4));
endfunction