function [ro,rou,rov,E] = inflow(ro,rou,rov,E, u)
ro(:)=u(1);
rou(:)=u(2);
rov(:)=u(3);
E(:)=u(4);
end

