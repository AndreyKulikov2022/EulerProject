function [ro0, rou0, rov0, E0]=InitializeLR(Nx, Ny, gamma, rol, roul, rovl, pl,ror, rour, rovr, pr )
El=energy(rol,roul,rovl,pl,gamma);
Er=energy(ror,rour,rovr,pr,gamma);
ro0=gpuArray([rol*ones(Ny+1,Nx/2+1),ror*ones(Ny+1,Nx/2)]);
rou0=gpuArray([roul*ones(Ny+1,Nx/2+1),rour*ones(Ny+1,Nx/2)]);
rov0=gpuArray([rovl*ones(Ny+1,Nx/2+1),rovr*ones(Ny+1,Nx/2)]);
E0=gpuArray([El*ones(Ny+1,Nx/2+1),Er*ones(Ny+1,Nx/2)]);