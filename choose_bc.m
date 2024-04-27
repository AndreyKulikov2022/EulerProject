function BC=choose_bc(bc)
if isnumeric(bc)
    BC=@(ro, rou, rov, E)inflow(ro, rou, rov, E,bc);
else
    switch bc
        case "zero_gradient"
            BC=@(ro, rou, rov, E)zero_gradient(ro, rou, rov, E);
        case "reflect_x"
            BC=@(ro, rou, rov, E)reflect_x(ro, rou, rov, E);
        case "reflect_y"
            BC=@(ro, rou, rov, E)reflect_y(ro, rou, rov, E);
    end
end
%Baically choose the correct value to append.