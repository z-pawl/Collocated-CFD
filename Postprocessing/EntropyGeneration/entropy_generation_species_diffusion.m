function entropy_generation = entropy_generation_species_diffusion(grid, energy_component, flow_component, species_manager_component)
    arguments
        grid (1,1) Grid2D
        energy_component (1,1) EnergyComponent
        flow_component (1,1) FlowComponent
        species_manager_component (1,1) SpeciesManagerComponent
    end
    
    %% Entropy generation due to species diffusion is calculated using the following equation
    % s_gen = -Jk * grad(μk/T); [W/(m^3*K)]
    % Where:
    % T - temperature [K], Jk - diffusive flux of species k [mol/(m^2*s)],
    % μk - chemical potential of species k [J/mol]

    % Chemical potential of species k is described using the following
    % formula
    % μk = μ0k + R * T * ln(pk / p0);
    % μk - chemical potential of species k [J/mol]
    % μ0k - standard chemical potential of species k [J/mol]
    % R - universal gas constant [J/(mol*K)]
    % T - temperature [K]
    % pi - partial pressure of species k [Pa]
    % p0 - reference pressure [Pa] = 1 bar = 10^5 Pa

    
    %% The following expression for the gradient of μk/T can be obtained
    % ∇(μk / T) = (T * d(μ0k)/dT - μ0k) / T^2 * ∇T + R * (∇P / P + ∇Yk / Yk + ∇Mmix / Mmix)
    % With additional symbols standing for:
    % P - pressure [Pa], Yk - mass fraction of species k [-], Mmix -
    % mixture molar mass [kg/mol]


    species_names = keys(species_manager_component.species);

    % Standard chemical potential is equal to the standard gibbs free
    % energy of formation as a polynomial of temperature [K]; [J/mol]
    standard_chemical_potentials_f = containers.Map;
    standard_chemical_potentials_f("H2") = [1.460925039083388e-12, -9.889856789569967e-09, 2.941472930603952e-05, -0.058611563335890, -1.049820977798776e+02, -2.928108543241390e+03];
    standard_chemical_potentials_f("CH4") = [4.941729056800309e-13, -3.259509426611343e-09, 1.343790053555073e-05, -0.061595600933205, -1.541029566151995e+02, -7.921666091023396e+04];
    standard_chemical_potentials_f("H2O") = [1.389558002141981e-12, -9.428229843354844e-09, 2.860942884155050e-05, -0.063782417707776, -1.602388475612041e+02, -2.451135763040180e+05];
    standard_chemical_potentials_f("CO") = [1.405590913550096e-12, -9.244343010357641e-09, 2.744105009595425e-05, -0.057491817887846, -1.719233246198568e+02, -1.135617562936047e+05];
    standard_chemical_potentials_f("CO2") = [6.639400309096650e-14, -3.083358883540760e-09, 1.892100949804213e-05, -0.066057266431854, -1.819515443924612e+02, -3.973209456773590e+05];

    % Product of the derivatives of the above functions w.r.t. temperature [K] and temperature [K]; [J/(mol)]
    T_times_std_chem_pot_f = containers.Map;
    T_times_std_chem_pot_f("H2") = conv([1 0], polyder(standard_chemical_potentials_f("H2")));
    T_times_std_chem_pot_f("CH4") = conv([1 0], polyder(standard_chemical_potentials_f("CH4")));
    T_times_std_chem_pot_f("H2O") = conv([1 0], polyder(standard_chemical_potentials_f("H2O")));
    T_times_std_chem_pot_f("CO") = conv([1 0], polyder(standard_chemical_potentials_f("CO")));
    T_times_std_chem_pot_f("CO2") = conv([1 0], polyder(standard_chemical_potentials_f("CO2")));


    % Gradient of following quantities have to be calculated: T -
    % temperature [K], P - pressure [Pa], Yi - mass fraction of species i [-],
    % Mmix - mixture molar mass [g/mol] (unit conversion isn't neccasary as
    % it will the gradient will be divided by the quantity)

    % Partial derivatives at cell centers (in particular direction) are 
    % calculated by taking the average of partial derivatives at faces 
    % normal to that direction

    % Coefficients for central differencing scheme
    [coeff_der_x, coeff_der_r] = central_differencing_scheme(grid);

    % Temperature gradient [K/m]
    [coeff_der_T_x, coeff_der_T_r] = energy_component.temp_bds.apply_boundary_condition_normal_derivative(coeff_der_x, coeff_der_r);
    [der_T_x, der_T_r] = evaluate_faces(coeff_der_T_x, coeff_der_T_r, energy_component.temp);
    gradient_T_x = (der_T_x(1:end-1,:) + der_T_x(2:end,:)) / 2;
    gradient_T_r = (der_T_r(:,1:end-1) + der_T_r(:,2:end)) / 2;
    clear coeff_der_T_x coeff_der_T_r der_T_x der_T_r;

    % Pressure gradient [Pa/m]
    [coeff_der_P_x, coeff_der_P_r] = flow_component.p_bds.apply_boundary_condition_normal_derivative(coeff_der_x, coeff_der_r);
    [der_P_x, der_P_r] = evaluate_faces(coeff_der_P_x, coeff_der_P_r, flow_component.p);
    gradient_P_x = (der_P_x(1:end-1,:) + der_P_x(2:end,:)) / 2;
    gradient_P_r = (der_P_r(:,1:end-1) + der_P_r(:,2:end)) / 2;
    clear coeff_der_P_x coeff_der_P_r der_P_x der_P_x;

    % Mass fractions gradients [1/m]
    gradient_Y_x = containers.Map;
    gradient_Y_r = containers.Map;
    for i = 1:numel(species_names)
        sp_name = species_names(i);
        sp = species_manager_component.species(sp_name);

        [coeff_der_Y_x, coeff_der_Y_R] = sp.species_bds.apply_boundary_condition_normal_derivative(coeff_der_x, coeff_der_r);
        [der_Y_x, der_Y_r] = evaluate_faces(coeff_der_Y_x, coeff_der_Y_R, sp.Y);
        gradient_Y_x(sp_name) = (der_Y_x(1:end-1,:) + der_Y_x(2:end,:)) / 2;
        gradient_Y_r(sp_name) = (der_Y_r(:,1:end-1) + der_Y_r(:,2:end)) / 2;
    end
    clear coeff_der_Y_x coeff_der_Y_R der_Y_x der_Y_r;

    % Mixture molar mass gradient [g/(mol*m)]
    % As the mixture molar mass is calculated using the following formula:
    % Mmix = (sum(Yi / Mi)) ^ (-1)
    % Its gradient can be calculated as follows:
    % ∇Mmix = -Mmix ^ 2 * sum(∇Yi / Mi)
    % Where: Yi - species i mass fraction [-], Mi - species i molar mass [g/mol]
    gradient_Mmix_x = zeros(grid.sz);
    gradient_Mmix_r = zeros(grid.sz);
    for i = 1:numel(species_names)
        sp_name = species_names(i);
        sp_M = species_manager_component.species(sp_name).M;
        gradient_Mmix_x = gradient_Mmix_x + gradient_Y_x(sp_name) / sp_M;
        gradient_Mmix_r = gradient_Mmix_r + gradient_Y_r(sp_name) / sp_M;
    end
    gradient_Mmix_x = -species_manager_component.M_mix .^ 2 .* gradient_Mmix_x;
    gradient_Mmix_r = -species_manager_component.M_mix .^ 2 .* gradient_Mmix_r;

    % Calculating ∇(μ/T) [J/(mol*K*m)]for each species
    gradient_chem_pot_over_T_x = containers.Map;
    gradient_chem_pot_over_T_r = containers.Map;
    for i = 1:numel(species_names)
        sp_name = species_names(i);
        sp = species_manager_component.species(sp_name);

        % The (T*dμ0/dT - μ0) / T^2 * ∇T term
        gradient_chem_pot_over_T_x_sp = (polyval(T_times_std_chem_pot_f(sp_name), energy_component.temp) - polyval(standard_chemical_potentials_f(sp_name), energy_component.temp)) ./ energy_component.temp .^ 2 .* gradient_T_x;
        gradient_chem_pot_over_T_r_sp = (polyval(T_times_std_chem_pot_f(sp_name), energy_component.temp) - polyval(standard_chemical_potentials_f(sp_name), energy_component.temp)) ./ energy_component.temp .^ 2 .* gradient_T_r;
    
        % The R * (∇P / P + ∇Yi / Yi + ∇Mmix / Mmix) term
        gradient_chem_pot_over_T_x_sp = gradient_chem_pot_over_T_x_sp + species_manager_component.R * ( ...
            gradient_P_x ./ flow_component.p ...
            + gradient_Y_x(sp_name) ./ (sp.Y) ...
            + gradient_Mmix_x ./ species_manager_component.M_mix);
        gradient_chem_pot_over_T_r_sp = gradient_chem_pot_over_T_r_sp + species_manager_component.R * ( ...
            gradient_P_r ./ flow_component.p ...
            + gradient_Y_r(sp_name) ./ (sp.Y) ...
            + gradient_Mmix_r ./ species_manager_component.M_mix);

        gradient_chem_pot_over_T_x(sp_name) = gradient_chem_pot_over_T_x_sp;
        gradient_chem_pot_over_T_r(sp_name) = gradient_chem_pot_over_T_r_sp;
    end

    %% Calculating the diffusion fluxes
    % The fluxes needed for the expression are molar fluxes, however first
    % the mass fluxes have to be calculated in order to apply mass
    % correction

    % Mass fluxes [kg/(m^2*s)] are given by the following expression:
    % J = -ϱ * D * ∇Y
    % Where: ϱ - density [kg/m^3], D - species diffusivity [m^2/s], Y -
    % species mass fraction [-]
    fluxes_x = containers.Map;
    fluxes_r = containers.Map;
    fluxes_sum_x = zeros(grid.sz);
    fluxes_sum_r = zeros(grid.sz);

    for i = 1:numel(species_names)
        sp_name = species_names(i);
        sp = species_manager_component.species(sp_name);

        flux_x_sp = -flow_component.rho .* sp.D .* gradient_Y_x(sp_name);
        flux_r_sp = -flow_component.rho .* sp.D .* gradient_Y_r(sp_name);

        fluxes_x(sp_name) = flux_x_sp;
        fluxes_r(sp_name) = flux_r_sp;
        fluxes_sum_x = fluxes_sum_x + flux_x_sp;
        fluxes_sum_r = fluxes_sum_r + flux_r_sp;
    end

    % Applying the correction of the form:
    % J* = J - Yi * sum(J)
    % And converting the mass fluxes to molar fluxes [mol/(m^2*s)]
    for i = 1:numel(species_names)
        sp_name = species_names(i);
        sp = species_manager_component.species(sp_name);

        M_sp = sp.M / 1000; % Species molar mass [kg/mol]

        fluxes_x(sp_name) = (fluxes_x(sp_name) - sp.Y .* fluxes_sum_x) / M_sp;
        fluxes_r(sp_name) = (fluxes_r(sp_name) - sp.Y .* fluxes_sum_r) / M_sp;
    end

    %% Calculating the entropy generation [W/(m^3*K)]
    entropy_generation = zeros(grid.sz);
    for i = 1:numel(species_names)
        sp_name = species_names(i);

        entropy_generation_sp = -(fluxes_x(sp_name) .* gradient_chem_pot_over_T_x(sp_name) + fluxes_r(sp_name) .* gradient_chem_pot_over_T_r(sp_name)); 
        % In cells where species mass fraction is smaller than 1e-8 the
        % entropy generation is assumed to be 0 in order to prevent wrong
        % results due to divisions by very small numbers
        mask = species_manager_component.species(sp_name).Y < 1e-8;
        entropy_generation_sp(mask) = 0;

        entropy_generation = entropy_generation + entropy_generation_sp;
    end
end