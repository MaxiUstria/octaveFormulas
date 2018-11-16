loglog(h,errorb2,'linewidth',2);
hold on;
loglog(h,errort2,'linewidth',2);
hold on;
loglog(h,errorh2,'linewidth',2);
hold on;
loglog(h,errorrk2,'linewidth',2);
hold on;
loglog(h,errore2,'linewidth',2);
hold on;
loglog(h,errordc2,'linewidth',2);
legend("Euler Hacia Atrás","Trapecio","Heun","RK-4","Euler Hacia Adelante","Diferencia Centrada");
hold on;
set (gca (), "xdir", "reverse");
xlabel("Paso","fontsize",16);
ylabel("Error Global","fontsize",16)'
