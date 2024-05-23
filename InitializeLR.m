function [ro0, rou0, rov0, E0]=InitializeLR(Nx, Ny, gamma, rol, roul, rovl, pl,ror, rour, rovr, pr,x0)
if nargin<12
x0=0.5;
end
El=energy(rol,roul,rovl,pl,gamma);
Er=energy(ror,rour,rovr,pr,gamma);
ro0=gpuArray([rol*ones(Ny+1,round(Nx*x0)+1),ror*ones(Ny+1,round(Nx*(1-x0)))]);
rou0=gpuArray([roul*ones(Ny+1,round(Nx*x0)+1),rour*ones(Ny+1,round(Nx*(1-x0)))]);
rov0=gpuArray([rovl*ones(Ny+1,round(Nx*x0)+1),rovr*ones(Ny+1,round(Nx*(1-x0)))]);
E0=gpuArray([El*ones(Ny+1,round(Nx*x0)+1),Er*ones(Ny+1,round(Nx*(1-x0)))]);