figure("Position",[0,0,2000,600])
hold on
top=polyshape([0,3,3,0],[0.5,0.5,0.8,0.8]+0.5);
bot=polyshape([0,0.6,0.6,3,3,0],[-0.5,-0.5,-0.3,-0.3,-0.8,-0.8]+0.5);
plot(union(top,bot),'FaceColor','black');
[X1,Y1]=meshgrid(0:0.05:3,[-0.25:0.05:0.45]+0.5);
X1=reshape(X1,[],1);
Y1=reshape(Y1,[],1);
[X2,Y2]=meshgrid(0:0.05:0.5,[-0.45:0.05:-0.25]+0.5);
X2=reshape(X2,[],1);
Y2=reshape(Y2,[],1);
X=[X1;X2];
Y=[Y1;Y2];
quiver(X,Y,0.04*ones(size(X,1),1),zeros(size(Y,1),1));
xlim([0,3]);
ylim([-0.8,0.8]+0.5);
grid minor;