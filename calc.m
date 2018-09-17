function [ x ] = calc(a,b,c,d)
  x(1,1) = c;
  x(1,2) = d;
  x(1,3) = ((-2*b)/(1 + a*a + b*b))*(c*d)
  x(1,4) = ((-2*a)/(1 + a*a + b*b))*(c*d)
endfunction