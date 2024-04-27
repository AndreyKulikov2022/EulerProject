function [Flux_x,Flux_y, BC_l, BC_b, BC_r, BC_t]=choose_method(method,bc_left, bc_bot, bc_right, bc_top)
%% Choose method and boundary conditions
switch method
    case "LF"
        Flux_x=@(ro, rou, rov, E, gamma, k, h)Euler_x_LF(ro, rou, rov, E, gamma, k, h);
        Flux_y=@(ro, rou, rov, E, gamma, k, h)Euler_y_LF(ro, rou, rov, E, gamma, k, h);
    case "HLLC"
        Flux_x=@(ro, rou, rov, E, gamma, k, h)Euler_x_HLLC(ro, rou, rov, E, gamma, k, h);
        Flux_y=@(ro, rou, rov, E, gamma, k, h)Euler_y_HLLC(ro, rou, rov, E, gamma, k, h);
end
BC_l=choose_bc(bc_left);
BC_b=choose_bc(bc_bot);
BC_r=choose_bc(bc_right);
BC_t=choose_bc(bc_top);
end

