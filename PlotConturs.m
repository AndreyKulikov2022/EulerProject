%% Density
figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
   contourf(X, Y, ro0, 'FaceColor',"white", 'LevelStep', 0.1)
%% U
figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
   contourf(X, Y, rou0./ro0, 'FaceColor',"white", 'LevelStep', 0.1)
%% Entropy
figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
contourf(X, Y, pressure(ro0,rou0,rov0,E0,gamma)./(ro0.^gamma), 'FaceColor',"white", 'LevelStep', 0.01)
%% Pressure
figure("Position",[0,0,1000*(Bx(end)-Bx(1))/(By(end)-By(1)),1000])
contourf(X, Y, pressure(ro0,rou0,rov0,E0,gamma), 'FaceColor',"white", 'LevelStep', 0.1)

%[X,Y]=meshgrid((Bx(1)+h/2):h:(Bx(2)-h/2),(By(1)+h/2):h:(By(2)-h/2));
%contourf(X, Y,[[ro_b(2:end-1,2:end-1),zeros(size(ro_b,1)-2,size(ro_r,2)-2)];[ro_l(2:end-1,2:end-1),ro_r(2:end-1,2:end-1)]], 'FaceColor',"white", 'LevelStep', 0.1);