%resolver ecuacion metodo octave
t= [0:0.001:2];
y0 = [-9 -8 3 2];
x= lsode ("fcn", y0, t);

