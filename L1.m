function l1 = L1(U,h)
%Find L2 norm of a function
l1=h*sum(abs(U),2);
end