%% Channel with a step
% Problem parameters
gamma = 1.4; % adiabatic gas constant
ro_i=1.4;
rou_i=3*ro_i;
rov_i=0;
p_i=1;
E_i=energy(ro_i,rou_i,rov_i,p_i,gamma);
Bx=[0,3]; % X boundaries ]
By=[0,1]; % Y boundaries
step_x=0.6; % step start position
step_y=0.2; % step height
T=4; % maximum time of computation
%% Method parameters
method="HLLC";
bc_x1="reflect_x"; % reflect from step
bc_x2="zero_gradient"; % zero gradient at the outlet
bc_bot="reflect_y";
bc_top="reflect_y";
[Flux_x, Flux_y, BC_x1, BC_b, BC_x2, BC_t]=choose_method(method,bc_x1, bc_bot, bc_x2, bc_top);
CFL=1/6;
h=0.0005; % space step
k=h*CFL;% time step
Nx=(Bx(2)-Bx(1))/h;
Ny=(By(2)-By(1))/h;
nx_step=step_x/h+1;
ny_step=step_y/h+1;
%[X,Y]=meshgrid(Bx(1):h:Bx(2),By(1):h:By(2));
%% Initialize data
ro0=gpuArray(ro_i*ones(Ny+1,Nx+1));
rou0=gpuArray(rou_i*ones(Ny+1,Nx+1));
rov0=gpuArray(rov_i*ones(Ny+1,Nx+1));
E0=gpuArray(E_i*ones(Ny+1,Nx+1));
% boundary helper arrays
ro_lb=ro_i*ones(ny_step-1,1);
rou_lb=rou_i*ones(ny_step-1,1);
rov_lb=rov_i*ones(ny_step-1,1);
E_lb=E_i*ones(ny_step-1,1);
ro_lt=ro_i*ones(Ny+2 - ny_step,1);
rou_lt=rou_i*ones(Ny+2 - ny_step,1);
rov_lt=rov_i*ones(Ny+2 - ny_step,1);
E_lt=E_i*ones(Ny+2 - ny_step,1);
%% Visualize
video2D=false;
if video2D
    figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
    hold on
    h_t=annotation('textbox',[.15 .6 .3 .3],'String','T = 0','FitBoxToText','on','FaceAlpha',0);
    im_ptr=imagesc('XData',Bx,'YData',By,'CData',rou0);
    xlim(Bx);
    ylim(By);
    clim([0.0,8]);
    colormap(turbo);
    %colormap(flipud(turbo));
    colorbar
%     testGIF='Project.gif';
%     F=getframe(gcf);
%     im=frame2im(F);
%     [imind,cm] = rgb2ind(im,128);
%     imwrite(imind,cm,testGIF,'gif','DelayTime',k, 'Loopcount',inf);
%     draw_count=3;
end
% Start calculation
t=0;
while t<T
    %% X - direction
    % Bottom part.
    [ro_rb, rou_rb, rov_rb, E_rb]=...
        BC_x1(ro0(1:ny_step-1,nx_step),rou0(1:ny_step-1,nx_step),rov0(1:ny_step-1,nx_step),E0(1:ny_step-1,nx_step));
    ro=[ro_lb,ro0(1:ny_step-1,1:nx_step),ro_rb];
    rou=[rou_lb,rou0(1:ny_step-1,1:nx_step),rou_rb];
    rov=[rov_lb,rov0(1:ny_step-1,1:nx_step),rov_rb];
    E=[E_lb,E0(1:ny_step-1,1:nx_step),E_rb];
    [ro0(1:ny_step-1,1:nx_step),rou0(1:ny_step-1,1:nx_step),rov0(1:ny_step-1,1:nx_step),E0(1:ny_step-1,1:nx_step)]=...
        solve_1d(ro, rou, rov, E, gamma, k, h, Flux_x,false);
    % Top part.
    [ro_rt, rou_rt, rov_rt, E_rt]=...
        BC_x2(ro0(ny_step:end,end),rou0(ny_step:end,end),rov0(ny_step:end,end),E0(ny_step:end,end));
    ro=[ro_lt,ro0(ny_step:end,1:end),ro_rt];
    rou=[rou_lt,rou0(ny_step:end,1:end),rou_rt];
    rov=[rov_lt,rov0(ny_step:end,1:end),rov_rt];
    E=[E_lt,E0(ny_step:end,1:end),E_rt];
    [ro0(ny_step:end,1:end),rou0(ny_step:end,1:end),rov0(ny_step:end,1:end),E0(ny_step:end,1:end)]=...
        solve_1d(ro, rou, rov, E, gamma, k, h, Flux_x,false);
    %% Y - direction
    % Left part.
    [ro_bb, rou_bb, rov_bb, E_bb]=BC_b(ro0(1,1:nx_step),rou0(1,1:nx_step),rov0(1,1:nx_step),E0(1,1:nx_step));
    [ro_tb, rou_tb, rov_tb, E_tb]=BC_t(ro0(end,1:nx_step),rou0(end,1:nx_step),rov0(end,1:nx_step),E0(end,1:nx_step));
    ro=[ro_bb;ro0(:,1:nx_step);ro_tb]';
    rou=[rou_bb;rou0(:,1:nx_step);rou_tb]';
    rov=[rov_bb;rov0(:,1:nx_step);rov_tb]';
    E=[E_bb;E0(:,1:nx_step);E_tb]';
    [ro0(:,1:nx_step),rou0(:,1:nx_step),rov0(:,1:nx_step),E0(:,1:nx_step)]=solve_1d(ro, rou, rov, E, gamma, k, h, Flux_y,1);
    % Right part.
    [ro_bb, rou_bb, rov_bb, E_bb]=BC_b(ro0(ny_step,nx_step+1:end),rou0(ny_step,nx_step+1:end),rov0(ny_step,nx_step+1:end),E0(ny_step,nx_step+1:end));
    [ro_tb, rou_tb, rov_tb, E_tb]=BC_t(ro0(end,nx_step+1:end),rou0(end,nx_step+1:end),rov0(end,nx_step+1:end),E0(end,nx_step+1:end));
    ro=[ro_bb;ro0(ny_step:end,nx_step+1:end);ro_tb]';
    rou=[rou_bb;rou0(ny_step:end,nx_step+1:end);rou_tb]';
    rov=[rov_bb;rov0(ny_step:end,nx_step+1:end);rov_tb]';
    E=[E_bb;E0(ny_step:end,nx_step+1:end);E_tb]';
    [ro0(ny_step:end,nx_step+1:end),rou0(ny_step:end,nx_step+1:end),rov0(ny_step:end,nx_step+1:end),E0(ny_step:end,nx_step+1:end)]=...
        solve_1d(ro, rou, rov, E, gamma, k, h, Flux_y,1);
    t=t+k;
    % Visualize.
    if video2D
        delete(im_ptr);
        h_t.String=['T = ',num2str(t)];
        im_ptr=imagesc('XData',Bx,'YData',By,'CData',rou0);
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
end

% figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
% hold on
% h_t=annotation('textbox',[.15 .6 .3 .3],'String','T = num2str(t)','FitBoxToText','on','FaceAlpha',0);
% im_ptr=imagesc('XData',Bx,'YData',By,'CData',pressure(ro0,rou0,rov0,E0,gamma)./(ro0.^gamma));
% xlim(Bx);
% ylim(By);
% clim([0.5,2]);
% colormap(turbo);
% %colormap(flipud(turbo));
% colorbar