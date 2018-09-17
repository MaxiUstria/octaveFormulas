clear;
%Superficie
x = [-11 :0.5:11];
y = x;
j = meshgrid(x,y);
k = j';
mesh(j,k,j.*k);
%Parametros
a = 0;
b = 2;
h = 0.000009;
y0 = [-9 -8 3 2];
%Euler hacia adelante
[t u v uP vP] = eu(a,b,y0,h);
hold on;
%Grafica la geodesica sobre la superficie
plot3(u,v,u.*v,"k",'linewidth',4);
%Diferencia Centrada
[t2 u2 v2 uP2 vP2] = difCen(a,b,y0,h);
hold on;
%Grafica de la geodesica sobre la superficie
plot3(u2,v2,(u2).*(v2),"r",'linewidth',4);
%resolver ecuacion metodo octave
p= lsode ("oct", y0, t);
hold on;
plot3(p(:,1),p(:,2),p(:,1).*p(:,2), 'b', "linewidth",8);