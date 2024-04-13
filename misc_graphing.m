distX=linspace(1,500,500)*h;
distY=linspace(1,201,201)*h;

f=figure(1);
 subplot(3,1,1)
 set(f,'position',[10,10,650,800])
 colormap 'jet'
 im1=imagesc(omega);
 axis tight manual equal
 set(0, 'defaulttextinterpreter','latex')
 title(strcat('Vorticity, t=',num2str(round(t,3))))
 cb=colorbar;
 caxis([-25,25])
 ylabel(cb, 'Vorticity (s^-1)')
 ylabel('Grid Point Number')
 xlabel('Grid Point Number')
 subplot(3,1,2)
 imagesc(temp)
 title(strcat('Temperature, t=',num2str(round(t,3))))
 cb=colorbar;
 caxis([300,400])
 ylabel(cb, "Temperature (C)")
 ylabel('Grid Point Number')
 xlabel('Grid Point Number')
 subplot(3,1,3)
 contourf(1:n,flip(1:m),psi,30)
 title(strcat('Stream Function, t=',num2str(round(t,3))))
 cb=colorbar;
 ylabel(cb, 'Stream Function (m^2/s)')
 ylabel('Grid Point Number')
 xlabel('Grid Point Number')
 

f2=figure(2);
set(f2,'position',[10,10,650,800])
subplot(3,1,1)
plot(distX,u(101,:))
xlabel("Distance from Left Border (m)")
ylabel('Horizontal Velocity u (m/s)')
title('Horizontal Velocity u along Horizontal Centerline, t=9.995')
xlim([distX(2),distX(499)]);
subplot(3,1,2)
plot(distX, v(101,:))
xlabel("Distance from Left Border (m)")
ylabel('Vertical Velocity v (m/s)')
title('Vertical Velocity v along Horizontal Centerline, t=9.995')
xlim([distX(2),distX(499)]);
subplot(3,1,3)
quiver(u,v)
xlim([200,250])
ylim([90,110])
title('Velocity Field On Small Area 0.15m from Right Edge of Golfball, t=10')
ylabel('Grid Point Number')
 xlabel('Grid Point Number')
f3=figure(3);
set(f3,'position',[10,10,650,800])
subplot(3,1,1)
plot(distY,u(:,100))
xlabel("Distance from Top Border (m)")
ylabel('Horizontal Velocity u (m/s)')
title('Horizontal Velocity u along Vertical Centerline of Golfball, t=9.995')
xlim([distY(2),distY(200)]);
subplot(3,1,2)
plot(distY, v(:,100))
xlabel("Distance from Left Border (m)")
ylabel('Vertical Velocity v (m/s)')
title('Vertical Velocity v along Vertical Centerline of Golfball, t=9.995')
xlim([distY(2),distY(200)]);
subplot(3,1,3)
quiver(u,v)
xlim([200,250])
ylim([90,110])
title('Velocity Field On Small Area 0.15m from Right Edge of Golfball, t=10')
ylabel('Grid Point Number')
 xlabel('Grid Point Number')
%}
time=linspace(1,10,645);
f4=figure(4);
set(f4,'position',[10,10,650,800])
subplot(3,1,1)
plot(distX,temp(101,:))
xlabel("Distance from Left Border (m)")
ylabel('Temperature (K)')
title('Temperature along Horizontal Centerline of Golfball, t=9.995')
xlim([distX(2),distX(499)]);
subplot(3,1,2)
plot(time,tempclose')
xlabel("Time (s)")
ylabel('Temperature (K)')
title('Temperature vs Time for a Point on the Horizontal Centerline and 1.46cm from Right Edge of Golfball')
subplot(3,1,3)
plot(time,tempfar')
xlabel("Time (s)")
ylabel('Temperature (K)')
title('Temperature vs Time for a Point on the Horizontal Centerline and 13.91cm from Right Edge of Golfball')

f5=figure(5);
plot(time,heat_boundary)
xlabel("Time (s)")
ylabel('Heat Flux (W/m)')
title('Heat Flux Exiting the Right Boundary over Time')
