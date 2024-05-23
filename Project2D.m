 %% Channel with a step (finite volume)
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
[Flux_x, Flux_y]=choose_method(method);
CFL=1/12;
h=0.0005; % space step
k=h*CFL;% time step
Nx=(Bx(2)-Bx(1))/h;
Ny=(By(2)-By(1))/h;
nx_step=round(step_x/h);
ny_step=round(step_y/h);
%[X,Y]=meshgrid(Bx(1):h:Bx(2),By(1):h:By(2));
%% Initialize data
%Bottomn left rectangle
ro_b=gpuArray((ro_i*ones(ny_step+2,nx_step+2)));
rou_b=gpuArray((rou_i*ones(ny_step+2,nx_step+2)));
rov_b=gpuArray((rov_i*ones(ny_step+2,nx_step+2)));
E_b=gpuArray((E_i*ones(ny_step+2,nx_step+2)));
%Top left rectangle
ro_l=gpuArray((ro_i*ones(Ny - ny_step+2,nx_step+2)));
rou_l=gpuArray((rou_i*ones(Ny - ny_step+2,nx_step+2)));
rov_l=gpuArray((rov_i*ones(Ny - ny_step+2,nx_step+2)));
E_l=gpuArray((E_i*ones(Ny - ny_step+2,nx_step+2)));
%Top right rectangle
ro_r=gpuArray((ro_i*ones(Ny - ny_step+2,Nx - nx_step+2)));
rou_r=gpuArray((rou_i*ones(Ny - ny_step+2,Nx - nx_step+2)));
rov_r=gpuArray((rov_i*ones(Ny - ny_step+2,Nx - nx_step+2)));
E_r=gpuArray((E_i*ones(Ny - ny_step+2,Nx - nx_step+2)));
%% Visualize
video2D=false;
if video2D
    figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
    hold on
    %[Xb,Yb]=meshgrid(Bx(1):h:step_x,By(1):h:step_y);
    h_t=annotation('textbox',[.15 .6 .3 .3],'String','T = 0','FitBoxToText','on','FaceAlpha',0);
    im_ptr_b=imagesc('XData',[Bx(1)+h/2,step_x-h/2],'YData',[By(1)+h/2,step_y-h/2],'CData',ro_b(2:end-1,2:end-1));
    im_ptr_l=imagesc('XData',[Bx(1)+h/2,step_x-h/2],'YData',[step_y+h/2, By(2)-h/2],'CData',ro_l(2:end-1,2:end-1));
    im_ptr_r=imagesc('XData',[step_x+h/2,Bx(2)-h/2],'YData',[step_y+h/2, By(2)-h/2],'CData',ro_r(2:end-1,2:end-1));
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
%     draw_count=8;
end
% Start calculation
t=0;
while t<T
    %% Boundary conditions
    %Bottomn left rectangle
    ro_b(end,2:end-1)=ro_l(2,2:end-1);
    rou_b(end,2:end-1)=rou_l(2,2:end-1);
    rov_b(end,2:end-1)=rov_l(2,2:end-1);
    E_b(end,2:end-1)=E_l(2,2:end-1);

    ro_b(1,2:end-1)=ro_b(2,2:end-1);
    rou_b(1,2:end-1)=rou_b(2,2:end-1);
    rov_b(1,2:end-1)=-rov_b(2,2:end-1);
    E_b(1,2:end-1)=E_b(2,2:end-1);

    ro_b(2:end-1,end)=ro_b(2:end-1,end-1);
    rou_b(2:end-1,end)=-rou_b(2:end-1,end-1);
    rov_b(2:end-1,end)=rov_b(2:end-1,end-1);
    E_b(2:end-1,end)=E_b(2:end-1,end-1);
    %Top left rectangle
    ro_l(1,2:end-1)=ro_b(end-1,2:end-1);
    rou_l(1,2:end-1)=rou_b(end-1,2:end-1);
    rov_l(1,2:end-1)=rov_b(end-1,2:end-1);
    E_l(1,2:end-1)=E_b(end-1,2:end-1);

    ro_l(end,2:end-1)=ro_l(end-1,2:end-1);
    rou_l(end,2:end-1)=rou_l(end-1,2:end-1);
    rov_l(end,2:end-1)=-rov_l(end-1,2:end-1);
    E_l(end,2:end-1)=E_l(end-1,2:end-1);

    ro_l(2:end-1,end)=ro_r(2:end-1,2);
    rou_l(2:end-1,end)=rou_r(2:end-1,2);
    rov_l(2:end-1,end)=rov_r(2:end-1,2);
    E_l(2:end-1,end)=E_r(2:end-1,2);
    %Top right rectangle
    ro_r(2:end-1,1)=ro_l(2:end-1,end-1);
    rou_r(2:end-1,1)=rou_l(2:end-1,end-1);
    rov_r(2:end-1,1)=rov_l(2:end-1,end-1);
    E_r(2:end-1,1)=E_l(2:end-1,end-1);

    ro_r(2:end-1,end)=ro_r(2:end-1,end-1);
    rou_r(2:end-1,end)=rou_r(2:end-1,end-1);
    rov_r(2:end-1,end)=rov_r(2:end-1,end-1);
    E_r(2:end-1,end)=E_r(2:end-1,end-1);

    ro_r(1,2:end-1)=ro_r(2,2:end-1);
    rou_r(1,2:end-1)=rou_r(2,2:end-1);
    rov_r(1,2:end-1)=-rov_r(2,2:end-1);
    E_r(1,2:end-1)=E_r(2,2:end-1);

    ro_r(end,2:end-1)=ro_r(end-1,2:end-1);
    rou_r(end,2:end-1)=rou_r(end-1,2:end-1);
    rov_r(end,2:end-1)=-rov_r(end-1,2:end-1);
    E_r(end,2:end-1)=E_r(end-1,2:end-1);
    %% Solve
    % Bottom part.
    [ro_b(2:end-1,2:end-1),rou_b(2:end-1,2:end-1),rov_b(2:end-1,2:end-1),E_b(2:end-1,2:end-1)]=...
        solve_2d(ro_b, rou_b, rov_b, E_b, gamma, k, h, Flux_x,Flux_y);
    % Top left part.
    [ro_l(2:end-1,2:end-1),rou_l(2:end-1,2:end-1),rov_l(2:end-1,2:end-1),E_l(2:end-1,2:end-1)]=...
        solve_2d(ro_l, rou_l, rov_l, E_l, gamma, k, h, Flux_x,Flux_y);
    % Top right part
    [ro_r(2:end-1,2:end-1),rou_r(2:end-1,2:end-1),rov_r(2:end-1,2:end-1),E_r(2:end-1,2:end-1)]=...
        solve_2d(ro_r, rou_r, rov_r, E_r, gamma, k, h, Flux_x,Flux_y);
    
    t=t+k;
    % Visualize.
    if video2D
        delete(im_ptr_b);
        delete(im_ptr_r);
        delete(im_ptr_l);
        h_t.String=['T = ',num2str(t)];
        im_ptr_b=imagesc('XData',[Bx(1)+h/2,step_x-h/2],'YData',[By(1)+h/2,step_y-h/2],'CData',ro_b(2:end-1,2:end-1));
        im_ptr_l=imagesc('XData',[Bx(1)+h/2,step_x-h/2],'YData',[step_y+h/2, By(2)-h/2],'CData',ro_l(2:end-1,2:end-1));
        im_ptr_r=imagesc('XData',[step_x+h/2,Bx(2)-h/2],'YData',[step_y+h/2, By(2)-h/2],'CData',ro_r(2:end-1,2:end-1));
        drawnow;
        %         if draw_count==8
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