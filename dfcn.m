function doct=df(x)
  doct(1,1) = 0;
  doct(2,1) = 0;
  doct(3,1) = (4*x(2)*x(1)/((1 + x(1)*x(1) + x(2)*x(2))^2))*(x(3)*x(4));
  doct(4,1) = (-2/((1 + x(1)*x(1) + x(2)*x(2))^2))*(x(3)*x(4))*(1 - x(1)*x(1) + x(2)*x(2));
  doct(1,2) = 0;
  doct(2,2) = 0;
  doct(3,2) = (-2/((1 + x(1)*x(1) + x(2)*x(2))^2))*(x(3)*x(4))*(1 + x(1)*x(1) - x(2)*x(2));;
  doct(4,2) = (4*x(2)*x(1)/((1 + x(1)*x(1) + x(2)*x(2))^2))*(x(3)*x(4));;
  doct(1,3) = 1;
  doct(2,3) = 0;
  doct(3,3) = (-2*x(2)/(1 + x(1)*x(1) + x(2)*x(2)))*(x(4));
  doct(4,3) = (-2*x(1)/(1 + x(1)*x(1) + x(2)*x(2)))*(x(4));
  doct(1,4) = 0;
  doct(2,4) = 1;
  doct(3,4) = (-2*x(2)/(1 + x(1)*x(1) + x(2)*x(2)))*(x(3));
  doct(4,4) = (-2*x(1)/(1 + x(1)*x(1) + x(2)*x(2)))*(x(3));
end