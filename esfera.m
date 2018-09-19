map= [0,0,0;0,0,0;0,0,0];
colormap(map);
r = 5;
fx = @(u,v) r.*cos (u) .* cos (v);
fy = @(u,v) r.*sin (u) .* cos (v);
fz = @(u,v) r.*sin (v);
ezmesh (fx, fy, fz, [0, 2*pi, -pi/2, pi/2], 30);
hold on
t = linspace(0,2.*pi,80);
cX = r.* cos(t);
cY = r.* sin(t);
plot3(cX,cY,'b','linewidth',3);
axis equal;
title (" ");
