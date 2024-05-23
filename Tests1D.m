%% Calculate 1D Riemann problem
% Problem parameters
gamma = 1.4; % adiabatic gas constant
Bx=[0,1]; % X boundaries
By=[0,1]; % Y boundaries
% Method parameters
direction="y";%1D direction
method="HLLC";
bc_left="zero_gradient";%[rol,roul,rovl,El];
bc_bot="zero_gradient";
bc_right="zero_gradient";
bc_top="zero_gradient";
[Flux_x, Flux_y, BC_l, BC_b, BC_r, BC_t]=choose_method(method,bc_left, bc_bot, bc_right, bc_top);
h=0.01; % space step
Nx=(Bx(2)-Bx(1))/h;
Ny=(By(2)-By(1))/h;
[X,Y]=meshgrid(Bx(1):h:Bx(2),By(1):h:By(2));
% Initialize values. %% Mind them being GPU arrays.
%rol=1;roul=0;rovl=0;pl=1;ror=0.125;rour=0;rovr=0;pr=0.1;T=0.2;x0=0.3;%Test1
rol=1;roul=-2;rovl=0;pl=0.4;ror=1;rour=2;rovr=0;pr=0.4;T=0.15;x0=0.5;%Test2
%rol=1;roul=0;rovl=0;pl=1000;ror=1;rour=0;rovr=0;pr=0.01;T=0.012;x0=0.5;%Test3
%rol=5.99924;roul=19.5975*rol;rovl=0;pl=460.894;ror=5.99242;rour=-6.19633*ror;rovr=0;pr=46.0950;T=0.035;x0=0.4;%Test4
%rol=1;roul=-19.59745;rovl=0;pl=1000;ror=1;rour=-19.59745;rovr=0;pr=0.01;T=0.012;x0=0.8;%Test5
if direction=="x"
[ro0, rou0, rov0, E0]=InitializeLR(Nx, Ny, gamma, rol, roul, 0, pl,ror, rour, 0, pr,x0);
elseif direction=="y"
[ro0, rou0, rov0, E0]=InitializeBT(Nx, Ny, gamma, rol, 0, roul, pl,ror, 0, rour, pr,x0);
end
%[ro0, rou0, rov0, E0]=InitializeCyl(Nx, Ny, 5/3, 1, 0, 0, 10,1, 0, 0, 0.1,X,Y);
%[ro0, rou0, rov0, E0]=InitializeBT(Nx, Ny, gamma, 1, 0, 0, 1, 0.125, 0, 0, 0.1);
Smax=max(max((abs(rou0./ro0)+sqrt(gamma*pressure(ro0,rou0,rov0,E0,gamma)./ro0)),[],"all"),max((abs(rov0./ro0)+sqrt(gamma*pressure(ro0,rou0,rov0,E0,gamma)./ro0)),[],"all"));
CFL=1/2;
k=h*CFL/Smax;% time step
% Visualization
show=true;
video2D=false;
video1Dx=false;
video1Dy=false;
if video2D
    figure("Position",[0,0,1500,1500])
    hold on
    im_ptr=imagesc('XData',Bx,'YData',By,'CData',ro0);
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
    %% X - direction
    % Apply boundary conditions using ghost points.
    [ro_lb, rou_lb, rov_lb, E_lb]=BC_l(ro0(:,1),rou0(:,1),rov0(:,1),E0(:,1));
    [ro_rb, rou_rb, rov_rb, E_rb]=BC_r(ro0(:,end),rou0(:,end),rov0(:,end),E0(:,end));
    ro=[ro_lb,ro0,ro_rb];
    rou=[rou_lb,rou0,rou_rb];
    rov=[rov_lb,rov0,rov_rb];
    E=[E_lb,E0,E_rb];
    % Obtain next step values.
    [ro0,rou0,rov0,E0]=solve_1d(ro, rou, rov, E, gamma, k, h, Flux_x,false);
    %% Y - direction
    % Apply boundary conditions using ghost points.
    [ro_bb, rou_bb, rov_bb, E_bb]=BC_b(ro0(1,:),rou0(1,:),rov0(1,:),E0(1,:));
    [ro_tb, rou_tb, rov_tb, E_tb]=BC_t(ro0(end,:),rou0(end,:),rov0(end,:),E0(end,:));
    % Obtain next step values.
    [ro0,rou0,rov0,E0]=solve_1d([ro_bb;ro0;ro_tb]', [rou_bb;rou0;rou_tb]', [rov_bb;rov0;rov_tb]', [E_bb;E0;E_tb]', gamma, k, h, Flux_y,1);

    t=t+k;
    if(t+k>T)
    k=T-t+0.000001;
    end
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
% [roa, ua,pa]=analit(1,0.125,1,0.1,gamma, X(1,:), T);
% hs(end+1)=h;
% L1s(end+1)=L1(reshape(ro0(1,:)-roa,1,[]),h);
% L2s(end+1)=L2(reshape(ro0(1,:)-roa,1,[]),h);
% Linfs(end+1)=Linf(reshape(ro0(1,:)-roa,1,[]));
if show && direction=="x"
    figure("Position",[0,0,500,500])
    hold on
    plot(X(10,:),ro0(1,:));
    %plot(X(1,:),roa);
    title("Density")
    figure("Position",[0,600,500,500])
    p=pressure(ro0, rou0, rov0, E0,gamma);
    plot(X(1,:),p(1,:));
    title("Pressure")
    figure("Position",[600,0,500,500])
    plot(X(1,:),rou0(1,:)./ro0(1,:));
    title("Velocity")
    figure("Position",[600,600,500,500])
    plot(X(1,:),p(1,:)/(gamma-1)./ro0(1,:));
    title("InternalEnergy")
end
if show && direction == "y"
 figure("Position",[1200,0,500,500])
    hold on
    plot(Y(:,1),ro0(:,1));
    %plot(X(1,:),roa);
    title("Density")
    figure("Position",[1200,600,500,500])
    p=pressure(ro0, rou0, rov0, E0,gamma);
    plot(Y(:,1),p(:,1));
    title("Pressure")
    figure("Position",[1800,0,500,500])
    plot(Y(:,1),rov0(:,1)./ro0(:,1));
    title("Velocity")
    figure("Position",[1800,600,500,500])
    plot(Y(:,1),p(:,1)/(gamma-1)./ro0(:,1));
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