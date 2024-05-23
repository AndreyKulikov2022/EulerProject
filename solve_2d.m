function [ro,rou,rov,E]=solve_2d(ro0, rou0, rov0, E0, gamma, k, h, Flux_x, Flux_y)
%% Calculates inner next step values in 1D
% ro0, rou0, rov0, E0 - initial values of the variables
% gamma - adiabatic gas constant
% k - time step
% h - space step
% Flux_x - Numeric x-flux function
% Flux_y - Numeric y-flux function

% %Calculate numeric fluxes
[F_ro_x,F_rou_x,F_rov_x,F_E_x]=Flux_x(ro0, rou0, rov0, E0, gamma, k, h);
[F_ro_y,F_rou_y,F_rov_y,F_E_y]=Flux_y(ro0', rou0', rov0', E0', gamma, k, h);
%Calculate variables at the next step
ro=ro0(2:end-1,2:end-1)-k/h*(F_ro_x(2:end-1,2:end)-F_ro_x(2:end-1,1:end-1))-k/h*(F_ro_y(2:end-1,2:end)'-F_ro_y(2:end-1,1:end-1)');
rou=rou0(2:end-1,2:end-1)-k/h*(F_rou_x(2:end-1,2:end)-F_rou_x(2:end-1,1:end-1))-k/h*(F_rou_y(2:end-1,2:end)'-F_rou_y(2:end-1,1:end-1)');
rov=rov0(2:end-1,2:end-1)-k/h*(F_rov_x(2:end-1,2:end)-F_rov_x(2:end-1,1:end-1))-k/h*(F_rov_y(2:end-1,2:end)'-F_rov_y(2:end-1,1:end-1)');
E=E0(2:end-1,2:end-1)-k/h*(F_E_x(2:end-1,2:end)-F_E_x(2:end-1,1:end-1))-k/h*(F_E_y(2:end-1,2:end)'-F_E_y(2:end-1,1:end-1)');
end