%% Compare results with accurate
% hs=[];
% L1s=[];
% L2s=[];
% Linfs=[];

h0=0.0005;
h=0.001;
ratio=h/h0;
sub_ind=gen_sub_ind(size(roa),ratio);

ro=reshape(roa(sub_ind),size(ro0,1),size(ro0,2));
rou=reshape(roua(sub_ind),size(ro0,1),size(ro0,2));
rov=reshape(rova(sub_ind),size(ro0,1),size(ro0,2));
E=reshape(Ea(sub_ind),size(ro0,1),size(ro0,2));

l1=h^2*(sum(abs(ro-ro0),"all")+sum(abs(rou-rou0),"all")+sum(abs(rov-rov0),"all")+sum(abs(E-E0),"all"));
l2=h*(sqrt(sum(abs(ro-ro0).^2,"all"))+sqrt(sum(abs(rou-rou0).^2,"all"))+sqrt(sum(abs(rov-rov0).^2,"all"))+sqrt(sum(abs(E-E0).^2,"all")));
linf=max([max(abs(ro-ro0),[],"all"),max(abs(rou-rou0),[],"all"),max(abs(rov-rov0),[],"all"),max(abs(E-E0),[],"all")]);

hs(end+1)=h;
L1s(end+1)=l1;
L2s(end+1)=l2;
Linfs(end+1)=linf;