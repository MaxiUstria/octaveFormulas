function p = promerror(u,v)
  i = 1;
  z = [];
  while i <= length(u)
    z(i) = sqrt(abs((u(i)*u(i))-(v(i)*v(i))));
    i = i+1;
  endwhile
  p = mean(z);
endfunction