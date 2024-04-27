function [ro0, rou0, rov0, E0]=InitializeCyl(Nx, Ny, gamma, roi, roui, rovi, pi,roo, rouo, rovo, po,X,Y)
Ei=energy(roi,roui,rovi,pi,gamma);
Eo=energy(roo,rouo,rovo,po,gamma);
ro0=gpuArray(roo*ones(Ny+1,Nx+1));
rou0=gpuArray(rouo*ones(Ny+1,Nx+1));
rov0=gpuArray(rovo*ones(Ny+1,Nx+1));
E0=gpuArray(Eo*ones(Ny+1,Nx+1));
% Change state in a circle in the middle
CXY=[(max(X,[],"all")+min(X,[],"all"))/2,(max(Y,[],"all")+min(Y,[],"all"))/2];
r=0.100001*min(max(X,[],"all")-min(X,[],"all"),max(Y,[],"all")-min(Y,[],"all"));
circle=polybuffer(CXY,"point",r);
id_in=zeros(size(X),"logical");
for i=1:size(X,2)
    id_in(:,i)=isinterior(circle,X(:,i),Y(:,i));
end
ro0(id_in)=roi;
rou0(id_in)=roui;
rov0(id_in)=rovi;
E0(id_in)=Ei;