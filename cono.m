r = 3;
c = 4*pi;
map= [0,0,0;0,0,0;0,0,0];
colormap(map);
fx = @(u,v) u.*cos (v) ;
fy = @(u,v) u.*sin (v) ;
fz = @(u,v) u - 4*pi;
ezmesh (fx, fy, fz, [0, 6*pi, 0,2*pi], 30);
hold on
plot3([0,0],[0,-6*pi],[-4*pi,2*pi],'g','linewidth',3);
hold on
j = linspace(0,2.*pi,80);
cX = c.* cos(j);
cY = c.* sin(j);
plot3(cX,cY,'r','linewidth',3);
hold on
t = linspace(0,6*pi,100);
cX = t.* cos(t);
cY = t.* sin(t);
cZ = t - 4*pi;
plot3(cX,cY,cZ,'b','linewidth',3);
axis equal;
title (" ");
