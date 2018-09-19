map= [0,0,0;0,0,0;0,0,0];
colormap(map);
r = pi;
c = pi;
fx = @(u,v) r.*cos (u) ;
fy = @(u,v) r.*sin (u) ;
fz = @(u,v) v;
ezmesh (fx, fy, fz, [0, 2*pi, -3*pi,3*pi], 30);
hold on
plot3([0,0],[-pi,-pi],[-3*pi,3*pi],'g','linewidth',3);
hold on
j = linspace(0,2.*pi,80);
cX = pi.* cos(j);
cY = pi.* sin(j);
plot3(cX,cY,'r','linewidth',3);
hold on
t = linspace(-3*pi,3*pi,50);
cX = r.* cos(t);
cY = r.* sin(t);
cZ = t;
plot3(cX,cY,cZ,'b','linewidth',2);
axis equal;
title (" ");
