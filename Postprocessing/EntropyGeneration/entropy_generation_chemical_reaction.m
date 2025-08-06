function entropy_generation = entropy_generation_chemical_reaction(grid, energy_component, flow_component, species_manager)
    arguments
        grid (1,1) Grid2D
        energy_component (1,1) EnergyComponent
        flow_component (1,1) FlowComponent
        species_manager (1,1) SpeciesManagerComponent
    end
    % Entropy generation due to a chemical reaction is calculated using the
    % following formula
    % s_gen = -(delta_G_r * r) / T; [W/(m^3*s*K)]
    % Where: T - temperature [K], r - reaction rate [mol/(m^3*s]
    % delta_G_r - Gibbs free energy change of a reaction [J/mol] which can
    % be calculated as following
    % delta_G_r = delta_G_r_0 + R * T * ln(Q)
    % delta_G_r_0 - standard Gibbs free energy change of a reaction [J/mol]
    % Q - reaction quotient [-]; R - universal gas constant = 8.314 [J/(mol*K)]

    % Reference pressure - 1 bar = 10^5 Pa
    p_ref = 1e5; % [Pa]
    
    % Mixture molar mass [g/mol]
    M_mix = species_manager.M_mix;

    % Involved species
    ch4 = species_manager.species("CH4");
    h2o = species_manager.species("H2O");
    co = species_manager.species("CO");
    co2 = species_manager.species("CO2");
    h2 = species_manager.species("H2");

    % A small number used to avoid division by zero
    eps = 1e-30;

    % Molar fractions [-]
    x_ch4 = ch4.Y .* M_mix / ch4.M + eps;
    x_h2o = h2o.Y .* M_mix / h2o.M + eps;
    x_co = co.Y .* M_mix / co.M + eps;
    x_co2 = co2.Y .* M_mix / co2.M + eps;
    x_h2 = h2.Y .* M_mix / h2.M + eps;


    % SMR reaction
    delta_G_0_SMR_f = [-5.726839566673235e-09, 3.112077001090684e-05, -0.064421568253950, -1.945144894332013e+02, 2.063669814811843e+05];
    delta_G_0_SMR = polyval(delta_G_0_SMR_f, energy_component.temp); % [J/mol]
    Q_SMR = (flow_component.p / p_ref) .^ 2 .* (x_co .* x_h2 .^ 3) ./ (x_ch4 .* x_h2o); % [-]
    delta_G_SMR = delta_G_0_SMR + species_manager.R * energy_component.temp .* log(Q_SMR);
    entropy_generation_smr = -delta_G_SMR .* species_manager.R_st ./ energy_component.temp;

    % WGS reaction
    delta_G_0_WGS = species_manager.delta_G; % [J/mol]
    Q_WGS = (x_co2 .* x_h2) ./ (x_co .* x_h2o); % [-]
    delta_G_WGS = delta_G_0_WGS + species_manager.R * energy_component.temp .* log(Q_WGS);
    entropy_generation_wgs = -delta_G_WGS .* species_manager.R_sh ./ energy_component.temp;

    entropy_generation = (entropy_generation_smr + entropy_generation_wgs);
end