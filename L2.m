function l2 = L2(U,h)
%Find L2 norm of a function
l2=sqrt(h*sum(U.^2,2));
end

