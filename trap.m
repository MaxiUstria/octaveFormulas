function [t,y]=trap(fcn,dfcn,y0,t0,tf,N)
  h=(tf-t0)/N;
  t=linspace(t0,tf,N+1);
  y=zeros(N+1,length(y0));
  y(1,:)=y0;
  for i=1:N
    x0=y(i,:)';
    x1=x0-inv(eye(length(y0))-(h/2)*feval(dfcn,x0))*(x0-(h/2)*(feval(fcn,x0)+feval(fcn,y(i,:)))'-y(i,:)');
    while (norm(x1-x0)>0.0001)
      x0=x1;
      x1=x0-inv(eye(length(y0))-(h/2)*feval(dfcn,x0))*(x0-(h/2)*(feval(fcn,x0)+feval(fcn,y(i,:)))'-y(i,:)');
    end
    y(i+1,:)=x1';
  end
end