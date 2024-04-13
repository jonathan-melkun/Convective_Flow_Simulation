% Initialize the m by n streamfunction array (psi)  
% of a specified cylinder radius (rad) centered at (circX)
close all;
uInf = 0.1;     %inflow velocity U_infinity
circX = 100;  %origin x-offset of cylinder 
rad = 25;     %cylinder radius
psi = zeros(201,500);
omega = zeros(201,500);
u = zeros(201,500);
v = zeros(201,500);
temp = zeros(201,500);
Re_grid=4;
nu=15.52e-6; %air
alpha=22.39e-6; %air
h=Re_grid*nu/uInf;
deltat=h/(4*uInf); 
[m,n] = size(psi);
k_air=26.24e-3;
tempclose=zeros(1,645); %1.49cm from right edge of cylinder
tempfar=zeros(1,645); %13.91cm from right edge of cylinder
heat_cylinder=zeros(1,645);
heat_boundary=zeros(1,645);
tempcount=1;
heatcount=1;

% Initialize psi array with parallel streamlines for the 'initial guess'
for i = 1:n
 psi(1:m,i) = uInf*h*flip((-m/2):(m/2-1));
end
% Set error tolerance and convergence flag to false;
ERR = 1e-3;
converge = false;
converge2 = false;


% Overwrite psi array with zeros at location of cylinder (surface and
% interior) and perform a 'Gauss-Seidel' interation until convergence flag
% is true (note the 'true/false' setting of the flag in the 'while loop'

%Initialize Inviscid Conditions
while ~converge
    converge = true;
 for i = 2:m-1
        for j = 2:n-1
            dist = sqrt((i-(m/2))^2 + (j-circX)^2);
            if (dist<=rad)
                psi(i,j) = 0;
            else
                resid = 1.5/4 * (psi(i-1,j) + psi(i+1,j) + psi(i,j-1) + psi(i, j+1) - 4*psi(i,j));
                if (resid/psi(i,j) > ERR)
                    converge = false;
                end
                psi(i,j) = psi(i,j) + resid;
            end
        end
 end
end
%Initialize Temperature Array
for i = 1:m
        for j = 1:n
            dist = sqrt((i-(m/2))^2 + (j-circX)^2);
            if (dist<=rad)
                temp(i,j)=400; %K
            else
                temp(i,j)=300;
            end
        end
end
% Initialize movie
prompt = 'Enter simulation movie filename: ';
    txt = input(prompt,"s");
    video = VideoWriter(txt);

video.FrameRate=30;
open(video);
% Main Time Step Loop
t=0;
tmax=10;
N=0;   
while t<tmax
    %Vorticity at the wall of cylinder + velocity + interior bulk vorticity
    omega_new=omega;
    temp_new=temp;
    for i = 2:m-1
        for j = 2:n-1
            dist = sqrt((i-(m/2))^2 + (j-circX)^2);
            if (dist<=rad)
                omega(i,j)=-2*((psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1))/h^2);
                omega_new(i,j)=-2*((psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1))/h^2);
                %temp
            else
                u(i,j)=((psi(i-1,j)-psi(i+1,j))/(2*h));
                v(i,j)=((psi(i,j-1)-psi(i,j+1))/(2*h));
            end
        end
    end
    % Advance Vorticity & Temp
    for i=2:m-1
        for j=2:n-1
            dist = sqrt((i-(m/2))^2 + (j-circX)^2);
            if (dist>rad)
                       
            if u(i,j)<0
                deltauw=u(i,j+1)*omega(i,j+1)-(u(i,j)*omega(i,j));
                deltaut=u(i,j)*(temp(i,j+1)-temp(i,j));
            else
                deltauw=u(i,j)*omega(i,j)-(u(i,j-1)*omega(i,j-1));
                deltaut=u(i,j)*(temp(i,j)-temp(i,j-1));
            end
            if v(i,j)<0
                deltavw=v(i-1,j)*omega(i-1,j)-(v(i,j)*omega(i,j));
                deltavt=v(i,j)*(temp(i-1,j)-temp(i,j));
            else
                deltavw=v(i,j)*omega(i,j)-(v(i+1,j)*omega(i+1,j));
                deltavt=v(i,j)*(temp(i,j)-temp(i+1,j));
            end
                omega_new(i,j) = omega(i,j) + deltat*(-deltauw/h-(deltavw/h)+nu*((omega(i+1,j)+omega(i-1,j)+omega(i,j+1)+omega(i,j-1)-4*omega(i,j))/h^2));
                temp_new(i,j) = temp(i,j) + deltat*(-deltaut/h-(deltavt/h)+alpha*((temp(i+1,j)+temp(i-1,j)+temp(i,j+1)+temp(i,j-1)-4*temp(i,j))/h^2));
            end
        end
    end
    omega=omega_new;
    temp=temp_new;
    % Calculate Bulk Vorticity to update psi array
    while ~converge2
        converge2=true;
    for i = m-1:-1:2
        for j = 2:n-1
            dist = sqrt((i-(m/2))^2 + (j-circX)^2);
            if (dist<=rad)
                psi(i,j)=0;
            else
            resid2 = 1.5/4 * (psi(i-1,j) + psi(i+1,j) + psi(i,j-1) + psi(i, j+1) + 4*(h^2)*omega(i,j) - 4*psi(i,j));
                if (resid2/psi(i,j) > ERR)
                    converge2 = false;
                end
                psi(i,j) = psi(i,j) + resid2;
            end
        end
    end
    
    end
    converge2=false;
    % Step 6: Calculate Outflow Boundary Condition
    for k=1:m
        psi(k,n)=2*psi(k,n-1)-psi(k,n-2);
        omega(k,n)=omega(k,n-1);
        temp(k,n)=2*temp(k,n-1)-temp(k,n-2);
    end
    
    %write the frame to the video
    % Plot data
    if (mod(N,5) == 0)
        f=figure(1);
        subplot(3,1,1)
        set(f,'position',[10,10,650,800])
        colormap 'jet'
        imagesc(omega)
        axis tight manual equal
        set(0, 'defaulttextinterpreter','latex')
        title(strcat('Vorticity, t=',num2str(round(t,3))))
        cb=colorbar;
        caxis([-25,25])
        ylabel(cb, 'Vorticity (s^-1)')
        subplot(3,1,2)
        imagesc(temp)
        title(strcat('Temperature, t=',num2str(round(t,3))))
        cb=colorbar;
        caxis([300,400])
        ylabel(cb, "Temperature (C)")
        subplot(3,1,3)
        contourf(1:n,flip(1:m),psi,30)
        title(strcat('Stream Function, t=',num2str(round(t,3))))
        cb=colorbar;
        ylabel(cb, 'Stream Function (m^2/s)')
        frame=getframe(f);
        writeVideo(video,frame);
    end
    N=N+1;
    t=deltat*N;
end
close(video)
 