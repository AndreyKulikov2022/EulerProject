function [F_ro,F_rou,F_rov,F_E]=Euler_y_LF(ro, rou, rov, E, gamma, k, h)
[f_ro,f_rou,f_rov,f_E]=Euler_y(ro,rou,rov,E,gamma);
F_ro=Flux_LF(ro,f_ro,h,k);
F_rou=Flux_LF(rou,f_rou,h,k);
F_rov=Flux_LF(rov,f_rov,h,k);
F_E=Flux_LF(E,f_E,h,k);
end