function [f_ro,f_rou,f_rov,f_E, p]=Euler_x(ro,rou,rov,E,gamma)
p=pressure(ro,rou,rov,E,gamma);
f_ro=rou;
f_rou=rou.^2./ro + p;
f_rov=rou.*rov./ro;
f_E=rou.*(E+p)./ro;
end