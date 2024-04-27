function p=pressure(ro,rou,rov,E,gamma)
p=(gamma-1)*(E-(rou.^2+rov.^2)./(2*ro));
end