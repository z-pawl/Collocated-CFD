% Class used for calculating the physical properties of gas mixture and reaction rates
classdef SpeciesManagerComponent < IComponent
    properties
        noi (1,1)                               % Number of iterations

        % Dependencies
        grid Grid2D {mustBeScalarOrEmpty}
        flow FlowComponent {mustBeScalarOrEmpty}
        energy EnergyComponent {mustBeScalarOrEmpty}
        species dictionary                      % Dictionary of chemical species

        Y_total (:,:) double                    % Sum of mass fractions of all species [-]
        M_mix (:,:) double                      % Molar mass of the mixture [g/mol]
        R_st (:,:) double                       % Rate of methane/steam reforming reaction [mol/(m^3*s)]
        R_sh (:,:) double                       % Rate of water-gas-shift reaction [mol/(m^3*s)]

        D_f (1,1) function_handle = @(x) 0      % Function used to calculate the effective mass diffusivity in the mixture for all species

        R double = 8.314472                     % Universal gas constant [J/(mol*K)]

        a double                                % Order of reaction with respect to methane [-]
        b double                                % Order of reaction with respect to water [-]
        A_st double                             % Arrhenius constant of the methane steam reforming reaction [mol/(g*s*Pa^(a+b))]
        E_a double                              % Activation energy of the methane steam reforming reaction [J/mol]
        delta_G double                          % Change of standard Gibbs free energy of water-gas-shift reaction [J/mol]
        w_cat double                            % Catalyst density [g/m^3]

        Jx_sum (:,:) double                     % Sum of diffusive fluxes in the x direction for all species [kg/s]  
        Jr_sum (:,:) double                     % Sum of diffusive fluxes in the r direction for all species [kg/s]

        % Solver settings
        relaxation_factor_R_st (1,1) double     % Relaxation factor of the steam-methane reforming reaction rate
        relaxation_factor_R_sh (1,1) double     % Relaxation factor of the water-gas-shift reaction rate
        inner_iters (1,1)                       % Number of inner iterations per outer iteration
        residual_history containers.Map         % Residual history of chemical species
    end
    methods
        % Constructor
        function obj = SpeciesManagerComponent(grid, flow, energy, D_function, a, b, A_st, E_a, delta_G, w_cat, relaxation_factor_R_st, relaxation_factor_R_sh, inner_iters)
            arguments
                grid (1,1) Grid2D
                flow (1,1) FlowComponent
                energy (1,1) EnergyComponent
                D_function (1,1) function_handle        % Function used to calculate species diffusivity in the mixture
                a (1,1) double                          % Rate of steam-methane reforming reaction with respect to methane [-]
                b (1,1) double                          % Rate of steam-methane reforming reaction with respect to water [-]
                A_st (1,1) double                       % Arrhenius constant [mol/(g*s*Pa^(a+b))]
                E_a (1,1) double                        % Activation energy [J/mol]
                delta_G (1,1) double                    % Change of standard Gibbs free energy of water-gas-shift reaction [J/mol]
                w_cat (1,1) double                      % Catalyst density [g/m^3]
                relaxation_factor_R_st (1,1) double     % Relaxation factor of the steam-methane reforming reaction rate
                relaxation_factor_R_sh (1,1) double     % Relaxation factor of the water-gas-shift reaction rate
                inner_iters (1,1) double                % Number of inner iterations per outer iteration
            end

            % Assigning the properties
            obj.noi = 0;

            obj.grid = grid;
            obj.flow = flow;
            obj.energy = energy;

            obj.Y_total = ones(grid.sz);
            obj.M_mix = ones(grid.sz);
            obj.R_st = zeros(grid.sz);
            obj.R_sh = zeros(grid.sz);
            
            obj.D_f = D_function;

            obj.a = a;
            obj.b = b;
            obj.A_st = A_st;
            obj.E_a = E_a;
            obj.delta_G = delta_G;
            obj.w_cat = w_cat;
            
            obj.Jx_sum = zeros(grid.sz + [1 0]);
            obj.Jr_sum = zeros(grid.sz + [0 1]);
            
            obj.relaxation_factor_R_st = relaxation_factor_R_st;
            obj.relaxation_factor_R_sh = relaxation_factor_R_sh;
            obj.inner_iters = inner_iters;
            obj.residual_history = containers.Map();
        end

        % Adds a species to the set of governed species
        function add_species(obj, name, species)
            arguments
                obj (1,1) SpeciesManagerComponent
                name (1,1) string
                species (1,1) Species
            end
            obj.species(name) = species;
            obj.residual_history(name) = zeros([max([obj.noi 10000]) 1]);
        end

        % Updates the physical properties and reaction rates, recalculates 
        % the sum of mass fractions for all species, molar mass of the 
        % mixture and the sum of diffusive fluxes of all species
        function update_properties(obj)
            species_components = values(obj.species);
            n_species = numel(species_components);

            % Sum of mass fractions [-]
            Y_total_new = zeros(obj.grid.sz);
            % Molar mass of the mixture [g/mol]
            % M_mix = 1 / sum(Yi/Mi); 
            % Yi - mass fraction of species i, Mi - molar mass of species i
            M_mix_new = zeros(obj.grid.sz);
            % Sum of diffusive fluxes of all species [kg/s]
            Jx_sum_new = zeros(obj.grid.sz + [1 0]);
            Jr_sum_new = zeros(obj.grid.sz + [0 1]);

            for i = 1:n_species
                % Clipping the negative mass fractions
                species_components(i).Y = clip(species_components(i).Y, 0, Inf);

                Y_total_new = Y_total_new + species_components(i).Y;
            end

            % Calculating the molar mass of the mixture and the sum of
            % diffusive fluxes
            for i = 1:n_species
                M_mix_new = M_mix_new + species_components(i).Y ./ (species_components(i).M  * Y_total_new);
                Jx_sum_new = Jx_sum_new + species_components(i).Jx;
                Jr_sum_new = Jr_sum_new + species_components(i).Jr;
            end
            M_mix_new = 1 ./ M_mix_new;

            % Calculating partial pressures for all species
            for i = 1:n_species
                species_components(i).p_partial = species_components(i).Y .* M_mix_new / species_components(i).M .* obj.flow.p; % [Pa]
            end

            % Assigning the calculated sum of mass fractions and the molar
            % mass of the mixture
            obj.Y_total = Y_total_new;
            obj.M_mix = M_mix_new;

            % % Updating the sum of diffusive fluxes of all species used by 
            % % the mass correction term for Fick's law based diffusion
            % obj.Jx_sum = Jx_sum_new;
            % obj.Jr_sum = Jr_sum_new;

            obj.update_species_diffusivity();
            obj.update_reaction_rates();
            obj.update_source_terms();
        end

        % Updates the effective mass diffusivity in the mixture for all species
        function update_species_diffusivity(obj)
            species_names = keys(obj.species);
            num_species = numel(species_names);

            D_species = obj.D_f();
            for i = 1:num_species
                name = species_names(i);
                obj.species(name).D = D_species(name); % [m^2/s]
            end
        end

        % Updates reaction rates of chemical reactions
        % SMR - rate is calculated using an empirical formula
        % WGS - rate is calculated assuming that the reaction is always in
        % equilibrium
        function update_reaction_rates(obj)            
            ch4_component = obj.species("CH4");
            h2o_component = obj.species("H2O");
            h2_component = obj.species("H2");
            co_component = obj.species("CO");
            co2_component = obj.species("CO2");

            % Correction: Y = Y* + b * v * R; b = volume * M / aP
            % Y - mass fraction [-], v - stoichiometric coefficient [-]
            % R - reaction rate [mol/(m^3*s)]
            % b [m^3*s/mol], volume [m^3], M - molar mass [kg/mol]
            % aP - central coefficient of the discretized transport equation of the species [kg/s]
            % Molar mass has to be divided by 1000 as its units are g/mol
            g_TO_kg = 1 / 1000;
            b_ch4 = obj.grid.volume * (ch4_component.M * g_TO_kg) ./ ch4_component.aP;
            b_h2o = obj.grid.volume * (h2o_component.M * g_TO_kg) ./ h2o_component.aP;
            b_h2 = obj.grid.volume * (h2_component.M * g_TO_kg) ./ h2_component.aP;
            b_co = obj.grid.volume * (co_component.M * g_TO_kg) ./ co_component.aP;
            b_co2 = obj.grid.volume * (co2_component.M * g_TO_kg) ./ co2_component.aP;

            % The methane steam reforming reaction rate is calculated explicitly
            % R_st = w_cat * A_st * exp(-E_a/(R*T)) * p_ch4^a * p_h2o^b [mol/(m^3*s)]
            % = w_cat * A_st * exp(-E_a/(R*T)) * p^(a+b) * X_ch4^a * X_h2o^b
            % = w_cat * A_st * exp(-E_a/(R*T)) * (p*M_mix)^(a+b) * Y_ch4^a * Y_h2o^b / (M_ch4^a*M_h2o^b)
            % w_cat - catalyst density [g/m^3], A_st - Arrhenius constant of the methane steam reforming reaction [mol/(g*s*Pa^(a+b))
            % E_a - activation energy of the methane steam reforming reaction [J/mol], R - universal gas constant [J/(mol*K)], T - temperature [K]
            % p_i - partial pressure of species i [Pa], p - pressure [Pa], X_i - molar fraction of species i [-], Y_i - mass fraction of species i [-]
            % M_mix - molar mass of the mixture [kg/mol], M_i - molar mass of species i [kg/mol] 
            G = obj.w_cat * obj.A_st .* exp(-obj.E_a ./ (obj.R .* obj.energy.temp)) .* (obj.flow.p .* obj.M_mix) .^ (obj.a + obj.b) / (ch4_component.M ^ obj.a * h2o_component.M ^ obj.b); % [mol/(m^3*s)]
            
            % Correction to the reaction rate is calculated and the rate is
            % updated using the underrelaxation factor
            R_st_corr = G .* ch4_component.Y .^ obj.a .* h2o_component.Y .^ obj.b ./ obj.Y_total .^ (obj.a + obj.b) - obj.R_st;
            obj.R_st = obj.R_st + obj.relaxation_factor_R_st * R_st_corr;
        

            % The water-gas-shift reaction rate is calculated using a rate
            % correction approach, the rate is relaxed to prevent oscillations

            % The equation for equilibrium constant is formulated using
            % mass fractions which are substituted by the correction formula
            % A quadratic equation is then obtained

            % The equilibrium constant [-] - it is multiplied by the molar
            % masses of species taking part in the reaction (units don't
            % have to be changed as they cancel one another)
            K_sh = exp(-obj.delta_G ./ (obj.R * obj.energy.temp)); 
            K_sh = K_sh * ((co2_component.M * h2_component.M) / (co_component.M * h2o_component.M)); % [-]

            % The coefficients of the quadratic equation
            A = K_sh .* b_co .* b_h2o - b_co2 .* b_h2; % [m^6*s^2/mol^2]
            B = -K_sh .* (b_co .* h2o_component.Y + b_h2o .* co_component.Y) - (b_co2 .* h2_component.Y + b_h2 .* co2_component.Y); % [m^3*s/mol]
            C = K_sh .* co_component.Y .*  h2o_component.Y - co2_component.Y .* h2_component.Y; % [-]

            B_2A = -B ./ (2 * A); % [mol/(m^3*s)]
            delta_2A = sqrt(B .* B - 4 * A .* C) ./ (2 * A); % [mol/(m^3*s)]

            % Updating the reaction rate using the relaxation factor
            R_sh_corr = B_2A - delta_2A; % [mol/(m^3*s)]
            obj.R_sh = obj.R_sh + obj.relaxation_factor_R_sh * R_sh_corr;
            % At least one species from each side of the reaction is
            % missing - the reaction rate is zero
            mask = (co_component.Y == 0 | h2o_component.Y == 0) & (co2_component.Y == 0 | h2_component.Y == 0);
            obj.R_sh(mask) = 0;
        end

        % Updates source terms for all species
        function update_source_terms(obj)
            sz = obj.grid.sz;

            % Stoichiometric coefficients for each species in each reaction [-]
            % Format: species_name: [ν_st, ν_sh] — stoichiometric coefficients
            % Reaction 1 (SMR): CH₄ + H₂O → CO + 3H₂
            % Reaction 2 (WGS): CO + H₂O ⇌ CO₂ + H₂

            stoich = {
              "H2", [3, 1];
              "CO", [1, -1];
              "CO2", [0 1];
              "CH4", [-1 0];
              "H2O", [-1 -1];
            };

            % The constant source term [kg/(m^3*s)]:
            % src_i = M_i * (v_st_i * R_st + v_sh_i * R_sh)
            % M_i - molar mass of species i [kg/mol] - values have to be
            % multiplied by 1/1000 as molar mass is in g/mol
            % v_st_i - stoichiometric coefficient of species i in the MSR reaction [-]; R_st - MSR reaction rate [mol/(m^3*s)]
            % v_sh_i - stoichiometric coefficient of species i in the WGS reaction [-]; R_sh - WGS reaction rate [mol/(m^3*s)]
            g_TO_kg = 1 / 1000;
            for i = 1:size(stoich,1)
                name = stoich{i, 1};
                coeffs = stoich{i, 2};

                sp = obj.species(name);
                sp.src = (sp.M * g_TO_kg) * (coeffs(1) * obj.R_st + coeffs(2) * obj.R_sh);
                % sp.src_lin = zeros(sz); - not changed anywhere else
            end
        end

        % Performs a single iteration and returns whether every species 
        % component has converged
        function converged = iterate(obj)
            converged = true;
            obj.noi = obj.noi + 1;

            species_components = values(obj.species);
            species_names = keys(obj.species);

            for i = 1:obj.inner_iters
                % Updating properties and performing an inner iteration for every species
                for j = 1:numel(species_components)
                    species_components(j).inner_iteration();
                end
                
                obj.update_properties();
            end

            % Updating properties
            obj.update_properties();

            % Calculating the residual and checking convergence
            for j = 1:numel(species_components)
                res = species_components(j).calculate_scaled_residual();
                res_hist = obj.residual_history(species_names(j));
                res_hist(obj.noi) = res;
                obj.residual_history(species_names(j)) = res_hist;
                converged = converged && res < species_components(j).tol;
            end
        end
    end
end