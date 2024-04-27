function E = energy(ro,rou,rov,p,gamma)
E=p/(gamma-1)+(rou.^2+rov.^2)./(2*ro);
end

