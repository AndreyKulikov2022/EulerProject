function [ro,rou,rov,E]=solve_1d(ro0, rou0, rov0, E0, gamma, k, h, Flux, transp_flag)
%% Calculates inner next step values in 1D
% ro0, rou0, rov0, E0 - initial values of the variables
% gamma - adiabatic gas constant
% k - time step
% h - space step
% flux - analytical flux function
% Flux - Numeric flux function

% %Calculate analyticall fluxes at each point
% [f_ro,f_rou,f_rov,f_E, p]=flux(ro0,rou0,rov0,E0,gamma);
% %Calculate numeric fluxes
[F_ro,F_rou,F_rov,F_E]=Flux(ro0, rou0, rov0, E0, gamma, k, h);
% F_ro=Flux(ro0,f_ro,p,h,k);
% F_rou=Flux(rou0,f_rou,p,h,k);
% F_rov=Flux(rov0,f_rov,p,h,k);
% F_E=Flux(E0,f_E,p,h,k);
%Calculate variables at the next step
ro=ro0(:,2:end-1)-k/h*(F_ro(:,2:end)-F_ro(:,1:end-1));
rou=rou0(:,2:end-1)-k/h*(F_rou(:,2:end)-F_rou(:,1:end-1));
rov=rov0(:,2:end-1)-k/h*(F_rov(:,2:end)-F_rov(:,1:end-1));
E=E0(:,2:end-1)-k/h*(F_E(:,2:end)-F_E(:,1:end-1));
if transp_flag
ro=ro';
rou=rou';
rov=rov';
E=E';
end
end

