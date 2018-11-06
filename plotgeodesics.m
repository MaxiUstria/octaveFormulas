xlabel('u(t)');
ylabel('v(t)');
scatter(u,v,50, "r");
hold on
scatter(u2,v2,50,"b");
hold on
plot(p(:,1),p(:,2),"k","linewidth",2);
pbaspect([1 1.3]);
hold on