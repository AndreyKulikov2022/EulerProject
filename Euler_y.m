function [f_ro,f_rou,f_rov,f_E, p]=Euler_y(ro,rou,rov,E,gamma)
p=pressure(ro,rou,rov,E,gamma);
f_ro=rov;
f_rou=rou.*rov./ro;
f_rov=rov.^2./ro + p;
f_E=rov.*(E+p)./ro;
end