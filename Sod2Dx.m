%% Calculate 1D Riemann problem
% Problem parameters
gamma = 1.4; % adiabatic gas constant
Bx=[0,1]; % X boundaries
By=[0,1]; % Y boundaries
T=0.2; % maximum time of computation
% Method parameters
method="HLLC";
[Flux_x, Flux_y]=choose_method(method);
CFL=1/4;
h=0.001; % space step
k=h*CFL;% time step
Nx=(Bx(2)-Bx(1))/h;
Ny=(By(2)-By(1))/h;
[X,Y]=meshgrid((Bx(1)+h/2):h:(Bx(2)-h/2),(By(1)+h/2):h:(By(2)-h/2));
rol=1;roul=0;rovl=0;pl=1;
ror=0.125;rour=0;rovr=0;pr=0.1;
El=energy(rol,roul,rovl,pl,gamma);
Er=energy(ror,rour,rovr,pr,gamma);
ro0=gpuArray([rol*ones(Ny+2,Nx/2+1),ror*ones(Ny+2,Nx/2+1)]);
rou0=gpuArray([roul*ones(Ny+2,Nx/2+1),rour*ones(Ny+2,Nx/2+1)]);
rov0=gpuArray([rovl*ones(Ny+2,Nx/2+1),rovr*ones(Ny+2,Nx/2+1)]);
E0=gpuArray([El*ones(Ny+2,Nx/2+1),Er*ones(Ny+2,Nx/2+1)]);
% Visualization
show=true;
video2D=false;
video1Dx=false;
video1Dy=false;
if video2D
    figure("Position",[0,0,1500,1500])
    hold on
    im_ptr=imagesc('XData',Bx,'YData',By,'CData',ro0(2:end-1,2:end-1));
    xlim(Bx);
    ylim(By);
    clim([0.0,5]);
    colormap(turbo);
    %colormap(flipud(turbo));
    colorbar
%     testGIF='MyVideo.gif';
%     F=getframe(gcf);
%     im=frame2im(F);
%     [imind,cm] = rgb2ind(im,128);
%     imwrite(imind,cm,testGIF,'gif','DelayTime',k, 'Loopcount',inf);
%     draw_count=3;
end

if video1Dx
    figure("Position",[0,0,500,500])
    hold on;
    xlim(Bx);
    ylim([min(ro0,[],"all")-0.001, max(ro0,[],"all")+0.001]);
    im1d_ptr=plot(X(1,:),ro0(1,:),'b');
    title("Density")
end

if video1Dy
    figure("Position",[0,0,500,500])
    hold on;
    xlim(Bx);
    ylim([min(ro0,[],"all")-0.001, max(ro0,[],"all")+0.001]);
    im1d_ptr=plot(Y(:,1),ro0(:,1),'b');
    title("Density")
end
% Start calculation
t=0;
while t<T
    %% BC
    ro0(2:end-1,1)=ro0(2:end-1,2);
    rou0(2:end-1,1)=rou0(2:end-1,2);
    rov0(2:end-1,1)=rov0(2:end-1,2);
    E0(2:end-1,1)=E0(2:end-1,2);

    ro0(2:end-1,end)=ro0(2:end-1,end-1);
    rou0(2:end-1,end)=-rou0(2:end-1,end-1);
    rov0(2:end-1,end)=rov0(2:end-1,end-1);
    E0(2:end-1,end)=E0(2:end-1,end-1);

    ro0(1,2:end-1)=ro0(2,2:end-1);
    rou0(1,2:end-1)=rou0(2,2:end-1);
    rov0(1,2:end-1)=rov0(2,2:end-1);
    E0(1,2:end-1)=E0(2,2:end-1);

    ro0(end,2:end-1)=ro0(end-1,2:end-1);
    rou0(end,2:end-1)=rou0(end-1,2:end-1);
    rov0(end,2:end-1)=rov0(end-1,2:end-1);
    E0(end,2:end-1)=E0(end-1,2:end-1);
   
    %% Solve
     [ro0(2:end-1,2:end-1),rou0(2:end-1,2:end-1),rov0(2:end-1,2:end-1),E0(2:end-1,2:end-1)]=...
        solve_2d(ro0, rou0, rov0, E0, gamma, k, h, Flux_x,Flux_y);
    t=t+k;
    % Visualize.
    if video2D
        delete(im_ptr);
        im_ptr=imagesc('XData',Bx,'YData',By,'CData',ro0);
        drawnow;
%         if draw_count==3
%         F=getframe(gcf);
%         im=frame2im(F);
%         [imind,cm] = rgb2ind(im,256);
%         imwrite(imind,cm,testGIF,'gif','DelayTime',k,'WriteMode','append');
%         draw_count=0;
%         else
%             draw_count=draw_count+1;
%         end
    end
    if video1Dx
        delete(im1d_ptr);
        im1d_ptr=plot(X(1,:),ro0(1,:),'b');
        drawnow;
    end
    if video1Dy
        delete(im1d_ptr);
        im1d_ptr=plot(Y(:,1),ro0(:,1),'b');
        drawnow;
    end
end
% 
[roa, ua,pa]=analit(1,0.125,1,0.1,gamma, X(1,:), T);
L1(reshape(ro0(2,2:end-1)-roa,1,[]),h)
if show
    figure("Position",[0,0,500,500])
    hold on
    plot(X(1,:),ro0(2,2:end-1));
    plot(X(1,:),roa);
    title("Density")
    figure("Position",[0,600,500,500])
    p=pressure(ro0, rou0, rov0, E0,gamma);
    plot(X(1,:),p(1,2:end-1));
    title("Pressure")
    figure("Position",[600,0,500,500])
    plot(X(1,:),rou0(1,2:end-1)./ro0(1,2:end-1));
    title("Velocity")
    figure("Position",[600,600,500,500])
    plot(X(1,:),p(1,2:end-1)/(gamma-1)./ro0(1,2:end-1));
    title("InternalEnergy")
end

function [ro, u,p]=analit(ro1,ro5,p1,p5,gamma, x, T)
ro=zeros(1,length(x));
u=zeros(1,length(x));
p=zeros(1,length(x));
c1=sqrt(gamma*p1/ro1);
c2=sqrt(gamma*p5/ro5);
x2=0.5-c1*T;
G=(gamma-1)/(gamma+1);
b=(gamma-1)/2/gamma;
p3=fsolve(@(p)(p1^b-p^b)*sqrt(((1-G^2)*p1^(1/gamma))/(G^2*p1))-...
    (p-p5)*sqrt((1-G)/(ro5*(p+G*p5))),(p1+p5)/2);
u3=(p1^b-p3^b)*sqrt(((1-G^2)*p1^(1/gamma))/(G^2*p1));
u4=(p3-p5)*sqrt((1-G)/(ro5*(p3+G*p5)));
p4=p3;
ro4=ro5*(p4+G*p5)/(p5+G*p4);
ro3=ro1*(p3/p1)^(1/gamma);
x3=0.5+(u3-sqrt(gamma*p3/ro3))*T;
x4=0.5+u3*T;
x5=0.5+c2*sqrt((gamma+1)*p4/2/gamma/p5+b)*T;
%formula for u2,p2,ro2
u2=2/(gamma+1)*(c1+(x(x>=x2&x<=x3)-0.5)/T);
ro2=ro1*(2/(gamma+1)-G/c1*(x(x>=x2&x<=x3)-0.5)/T).^(2/(gamma-1));%ro1*(1-(gamma-1)/2*u2/c1).^(2/(gamma-1));
p2=p1*(1-(gamma-1)/2*u2/c1).^(2*gamma/(gamma-1));

ro(x<x2)=ro1;
ro(x>=x2&x<=x3)=ro2;
ro(x>x3)=ro3;
ro(x>x4)=ro4;
ro(x>x5)=ro5;

u(x>=x2&x<=x3)=u2;
u(x>x3)=u3;
u(x>x4)=u4;
u(x>x5)=0;

p(x<x2)=p1;
p(x>=x2&x<=x3)=p2;
p(x>x3)=p3;
p(x>x4)=p4;
p(x>x5)=p5;

end