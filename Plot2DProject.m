%% Density
figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
    hold on
    h_t=annotation('textbox',[.15 .6 .3 .3],'String','T = 4','FitBoxToText','on','FaceAlpha',0);
    im_ptr=imagesc('XData',Bx,'YData',By,'CData',ro0);
    plot(polyshape([step_x,Bx(end),Bx(end),step_x],[By(1),By(1),step_y,step_y]),'FaceColor','black')
    xlim(Bx);
    ylim(By);
    clim([0.0,8]);
    colormap(turbo);
    %colormap(flipud(turbo));
    colorbar
%% U
figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
    hold on
    h_t=annotation('textbox',[.15 .6 .3 .3],'String','T = 4','FitBoxToText','on','FaceAlpha',0);
    im_ptr=imagesc('XData',Bx,'YData',By,'CData',rou0./ro0);
    plot(polyshape([step_x,Bx(end),Bx(end),step_x],[By(1),By(1),step_y,step_y]),'FaceColor','black')
    xlim(Bx);
    ylim(By);
    clim([0.0,3.5]);
    colormap(turbo);
    %colormap(flipud(turbo));
    colorbar
%% V
figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
    hold on
    h_t=annotation('textbox',[.15 .6 .3 .3],'String','T = 4','FitBoxToText','on','FaceAlpha',0);
    im_ptr=imagesc('XData',Bx,'YData',By,'CData',rov0./ro0);
    xlim(Bx);
    ylim(By);
plot(polyshape([step_x,Bx(end),Bx(end),step_x],[By(1),By(1),step_y,step_y]),'FaceColor','black')
    clim([0.0,1.5]);
    colormap(turbo);
    %colormap(flipud(turbo));
    colorbar
%% Entropy
figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
hold on
h_t=annotation('textbox',[.15 .6 .3 .3],'String','T = 4','FitBoxToText','on','FaceAlpha',0);
im_ptr=imagesc('XData',Bx,'YData',By,'CData',pressure(ro0,rou0,rov0,E0,gamma)./(ro0.^gamma));
plot(polyshape([step_x,Bx(end),Bx(end),step_x],[By(1),By(1),step_y,step_y]),'FaceColor','black')
xlim(Bx);
ylim(By);
clim([0.5,2]);
colormap(turbo);
%colormap(flipud(turbo));
colorbar
%% Pressure
figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
hold on
h_t=annotation('textbox',[.15 .6 .3 .3],'String','T = 4','FitBoxToText','on','FaceAlpha',0);
im_ptr=imagesc('XData',Bx,'YData',By,'CData',pressure(ro0,rou0,rov0,E0,gamma));
plot(polyshape([step_x,Bx(end),Bx(end),step_x],[By(1),By(1),step_y,step_y]),'FaceColor','black')
xlim(Bx);
ylim(By);
clim([0.5,12]);
colormap(turbo);
%colormap(flipud(turbo));
colorbar
plot(polyshape([step_x,Bx(end),Bx(end),step_x],[By(1),By(1),step_y,step_y]),'FaceColor','black')
plot(polyshape([step_x,Bx(end),Bx(end),step_x],[By(1),By(1),step_y,step_y]),'FaceColor','black')
plot(polyshape([step_x,Bx(end),Bx(end),step_x],[By(1),By(1),step_y,step_y]),'FaceColor','black')
plot(polyshape([step_x,Bx(end),Bx(end),step_x],[By(1),By(1),step_y,step_y]),'FaceColor','black')