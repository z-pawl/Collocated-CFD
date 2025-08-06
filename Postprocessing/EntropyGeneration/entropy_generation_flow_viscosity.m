function entropy_generation = entropy_generation_flow_viscosity(grid, energy_component, flow_component)
    arguments
        grid (1,1) Grid2D
        energy_component (1,1) EnergyComponent
        flow_component (1,1) FlowComponent
    end
    % Entropy generation due to flow viscosity in a porous media is
    % described by the following equation:
    % s_gen = -(v * grad(p)) / T; [W/(m^3*s*K)]
    % Where: T - temperature [K], v - velocity [m/s]
    % p - pressure [Pa]

    % Calculation of the pressure gradient
    % Pressure at faces
    [coeff_p_x, coeff_p_r] = linear_interpolation_scheme(grid);
    [coeff_p_x, coeff_p_r] = flow_component.p_bds.apply_boundary_condition_value(coeff_p_x, coeff_p_r);
    [p_x, p_r] = evaluate_faces(coeff_p_x, coeff_p_r, flow_component.p); % [Pa]

    % Pressure gradient [Pa/m]
    grad_p_x = (p_x(2:end,:) - p_x(1:end-1,:)) ./ grid.dx;
    grad_p_r = (p_r(:,2:end) - p_r(:,1:end-1)) ./ grid.dr;

    entropy_generation = -(flow_component.vx .* grad_p_x + flow_component.vr .* grad_p_r) ./ energy_component.temp;
end