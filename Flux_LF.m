function F=Flux_LF(u,fu,h,k)
F=1/2*(fu(:,2:end)+fu(:,1:end-1))-(h/(2*k))*(u(:,2:end)-u(:,1:end-1));
end