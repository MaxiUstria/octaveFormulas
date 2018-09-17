function p = promerror(u,v)
  i = 1;
  z = [];
  while i <= length(u)
    z(1,i) = sqrt((u(1,i)^2)*(v(1,i)^2))
    i = i+1;
  endwhile
  p = mean(z);
endfunction