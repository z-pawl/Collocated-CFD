function entropy_generation = entropy_generation_flow_viscosity(grid, energy_component, flow_component, permeability)
    arguments
        grid (1,1) Grid2D
        energy_component (1,1) EnergyComponent
        flow_component (1,1) FlowComponent
        permeability (1,1) double
    end
    % Entropy generation due to flow viscosity in a porous media is
    % described by the following equation:
    % s_gen = visc * v^2 / (Kp * T); [W/(m^3*s*K)]
    % Where: T - temperature [K], v - velocity [m/s]
    % visc - dynamic viscosity [Pa*s], Kp - permeability [m^2]

    entropy_generation = flow_component.visc .* (flow_component.vx .^ 2 + flow_component.vr .^ 2) ./ (permeability * energy_component.temp);
end