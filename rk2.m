function [t, y] = rk2(fcn,t0,tf,y0,N)
h = (tf - t0) / N;
y(1,:) = y0;
t(1) = t0;
for i = 1 : N
    t(i+1) = t(i) + h;
    K1 = feval(fcn, y(i,:));
    K2 = (feval(fcn, y(i,:) + K1./2));
    y(i+1,:) = y(i,:) +  (K1+K2)/2;
end;