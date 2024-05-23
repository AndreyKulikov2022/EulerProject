function [ro0, rou0, rov0, E0]=InitializeBT(Nx, Ny, gamma, rob, roub, rovb, pb,rot, rout, rovt, pt, y0)
if nargin<12
y0=0.5;
end
Eb=energy(rob,roub,rovb,pb,gamma);
Et=energy(rot,rout,rovt,pt,gamma);
ro0=gpuArray([rob*ones(round(Ny*y0)+1,Nx+1);rot*ones(round(Ny*(1-y0)),Nx+1)]);
rou0=gpuArray([roub*ones(round(Ny*y0)+1,Nx+1);rout*ones(round(Ny*(1-y0)),Nx+1)]);
rov0=gpuArray([rovb*ones(round(Ny*y0)+1,Nx+1);rovt*ones(round(Ny*(1-y0)),Nx+1)]);
E0=gpuArray([Eb*ones(round(Ny*y0)+1,Nx+1);Et*ones(round(Ny*(1-y0)),Nx+1)]);