function [ro0, rou0, rov0, E0]=InitializeBT(Nx, Ny, gamma, rob, roub, rovb, pb,rot, rout, rovt, pt)
Eb=energy(rob,roub,rovb,pb,gamma);
Et=energy(rot,rout,rovt,pt,gamma);
ro0=gpuArray([rob*ones(Ny/2+1,Nx+1);rot*ones(Ny/2,Nx+1)]);
rou0=gpuArray([roub*ones(Ny/2+1,Nx+1);rout*ones(Ny/2,Nx+1)]);
rov0=gpuArray([rovb*ones(Ny/2+1,Nx+1);rovt*ones(Ny/2,Nx+1)]);
E0=gpuArray([Eb*ones(Ny/2+1,Nx+1);Et*ones(Ny/2,Nx+1)]);