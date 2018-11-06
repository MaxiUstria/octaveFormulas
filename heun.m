function [t, y] = heun(fcn,t0,tf,y0,N)
h = (tf - t0) / N;
y(1,:) = y0;
t(1) = t0;
for i = 1 :  N
    t(i+1) = t(i) + h;
    z = y(i,:) + h * feval(fcn,y(i,:));
    y(i+1,:) = y(i,:) + (h/2) * ( feval(fcn,y(i,:)) + feval(fcn,z) );
end;