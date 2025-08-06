classdef ThermophysicalProperties < handle
    properties
        % Dependencies
        grid Grid2D {mustBeScalarOrEmpty}
        flow FlowComponent {mustBeScalarOrEmpty}
        energy EnergyComponent {mustBeScalarOrEmpty}
        species_manager SpeciesManagerComponent {mustBeScalarOrEmpty}

        D_AB double             % Binary diffusion coefficients (without the part depending on temperature and pressure) []
        phi_visc double         % Coefficients used to calculate the mixture viscosity [-]

        solid_phase_thermal_conductivity (1,1) double % Solid phase thermal conductivity [W/(m*K)]
    end
    methods
        function obj = ThermophysicalProperties(grid, flow, energy, species_manager, solid_phase_thermal_conductivity)
            arguments
                grid (1,1) Grid2D
                flow (1,1) FlowComponent
                energy (1,1) EnergyComponent
                species_manager (1,1) SpeciesManagerComponent
                solid_phase_thermal_conductivity (1,1) double % Solid phase thermal conductivity [W/(m*K)]
            end

            obj.grid = grid;
            obj.flow = flow;
            obj.energy = energy;
            obj.species_manager = species_manager;
            obj.solid_phase_thermal_conductivity = solid_phase_thermal_conductivity;

            obj.update_values();
        end
        
        % Function updating the values of D_AB and phi_visc
        function update_values(obj)
            obj.calculate_D_AB();
            obj.calculate_phi_visc();
        end

        % Function calculating binary diffusivities (without the p and T
        % part) used to calculate the effective mass diffusivity in the mixture
        function calculate_D_AB(obj)
            C = 0.00143;  % Empirical constant
            species_list = obj.species_manager.species;
            species_list = values(species_list);
            num_species = numel(species_list);
        
            % Preallocate diffusivity matrix
            D_AB_temp = zeros(num_species);
        
            % Precompute properties
            M = zeros(num_species, 1);      % Molar masses [g/mol]
            V_D = zeros(num_species, 1);    % Diffusion volumes [?]
        
            for i = 1:num_species
                M(i) = species_list(i).M;
                V_D(i) = species_list(i).V_D;
            end
        
            % Compute binary diffusion coefficients
            for i = 1:num_species
                for j = 1:num_species
                    if i == j
                        D_AB_temp(i, j) = Inf;
                    else
                        M_AB = 2 / (1 / M(i) + 1 / M(j)); % [-]
                        sigma = (V_D(i)^(1/3) + V_D(j)^(1/3))^2; % [?]
                        D_AB_temp(i, j) = C / (sqrt(M_AB) * sigma);
                    end
                end
            end
        
            obj.D_AB = D_AB_temp;
        end

        % Function used to calculate phi used in the mixture dynamic
        % viscosity formula. Phi is calculated using the method of Herning and Zipperer [-]
        function calculate_phi_visc(obj)
            species_list = values(obj.species_manager.species); 
            num_species = numel(species_list);

            phi = zeros(num_species);

            for i = 1:num_species
                for j = 1:num_species
                    phi(i,j) = (species_list(j).M / species_list(i).M) ^ (0.5);
                end
            end

            obj.phi_visc = phi;
        end

        % Function used to calculate the effective mass diffusivity in the
        % mixture for all species [m^2/s]
        function D = D_function(obj)
            % Returns mixture-averaged diffusivities for all species
            % using the Wilke formula.

            species_struct = obj.species_manager.species;
            species_names = keys(species_struct);
            species_list = values(species_struct);
            num_species = numel(species_list);
        
            % Constants for unit conversion
            PA_TO_BAR = 1e-5;
            CM2S_TO_M2S = 1e-4;
        
            % Precompute molar fractions
            M_mix = obj.species_manager.M_mix; % Mixture molar mass [g/mol]
            Y_tot = obj.species_manager.Y_total; % Sum of mass fractions [-]
            molar_fractions = cell(num_species,1); % Molar fractions [-]
            
            for i = 1:num_species
                sp = species_list(i);
                molar_fractions{i} = (sp.Y ./ Y_tot) .* M_mix ./ sp.M; % [-]
            end
        
            % Allocate output map for diffusivities
            D = containers.Map();
        
            % Common parameters
            T = obj.energy.temp;     % Temperature field [K]
            P = obj.flow.p;          % Pressure field [Pa]
        
            % Compute diffusivities for each species
            for i = 1:num_species
                name_i = species_names(i);
                D_mix_inv = zeros(size(T));
        
                for j = 1:num_species
                    if j ~= i
                        X_j = molar_fractions{j};
                        D_ij = obj.D_AB(i, j);  % Binary diffusivity
                        D_mix_inv = D_mix_inv + X_j / D_ij;
                    end
                end
        
                D_i = (T .^ 1.75) ./ (P * PA_TO_BAR) ./ D_mix_inv;
                D_i = D_i * CM2S_TO_M2S;  % Convert to m²/s
                D_i = (1 - sqrt(1 - obj.flow.porosity)) .* D_i; % Calculating the effective diffusivity [m^2/s]
                D(name_i) = D_i;
            end
        end

        % Function used to calculate the density using the ideal gas law [kg/m^3]
        function rho = rho_function(obj)
            % The density is calculated using the ideal gas law
            % rho = p * M_mix / (R * T) [kg/m^3]
            % p - pressure [Pa], M_mix - mixture molar mass [kg/mol] - has to be divided by 1000 as the unit of molar mass is g/mol
            % R - universal gas constant [J/(mol*K)], T - temperature [K]
            g_TO_kg = 1 / 1000;
            rho = obj.flow.p .* (obj.species_manager.M_mix * g_TO_kg) ./ (obj.species_manager.R * obj.energy.temp);
        end
    
        % Function used to calculate the heat capacity [J/(kg*K)]
        function cp = cp_function(obj)
            species_map = obj.species_manager.species;
            species_list = values(species_map);
            num_species = numel(species_list);
            
            T = obj.energy.temp;  % Temperature field [K]
            M_mix = obj.species_manager.M_mix; % [g/mol]
            Y_tot = obj.species_manager.Y_total; % Sum of mass fractions [-]
        
            % Preallocate arrays
            cp = zeros(size(T));
            
            for i = 1:num_species
                sp = species_list(i);
                
                % Molar fraction: Y * M_mix / M_i, [-]
                X_i = (sp.Y ./ Y_tot) .* M_mix ./ sp.M;
                
                % Molar heat capacity [J/(mol*K)]
                cp_i = polyval(sp.Cp, T);
                
                cp = cp + X_i .* cp_i;
            end
        
            % Convert from molar to specific heat capacity
            % Unit of molar mass is g/mol, thus it has to be changed
            g_TO_kg = 1 / 1000;
            cp = cp .* (M_mix * g_TO_kg);
        end
    
        % Function used to calculate the thermal conductivity [W/(m*K)]
        function k = k_function(obj)
            % Nested helper function for phi [-]
            function phi = phi_f(Mi, Mj, visc_i, visc_j) 
                M_ratio = Mi / Mj; % [-]
                visc_ratio = visc_i ./ visc_j; % [-]

                phi = (1 + sqrt(visc_ratio) * M_ratio ^ (-0.25)) .^ 2 / sqrt(8 * (1 + M_ratio));
            end

            species_map = obj.species_manager.species;
            species_list = values(species_map);
            num_species = numel(species_list);

            T = obj.energy.temp;  % Temperature field [K]
            M_mix = obj.species_manager.M_mix; % Mixture molar mass [g/mol]
            Y_tot = obj.species_manager.Y_total; % Sum of mass fractions [-]
            sz = size(T);

            % Preallocating the arrays
            k = zeros(sz); % Mixture thermal conductivity [W/(m*K)]

            X_i = cell(num_species,1); % Molar fractions [-]
            visc_i = cell(num_species, 1); % Dynamic viscosities of pure substances [Pa*s]
            k_i = cell(num_species, 1); % Thermal conductivities of pure substances [W/(m*K)]

            for i = 1:num_species
                sp = species_list(i);
                
                % Molar fraction: Y * M_mix / M_i [-]
                X_i{i} = (sp.Y ./ Y_tot) .* M_mix ./ sp.M;
                
                % Dynamic viscosity [Pa*s]
                visc_i{i} = polyval(sp.visc, T);

                % Thermal conductivity [W/(m*K)]
                k_i{i} = polyval(sp.k, T);
            end

            for i = 1:num_species
                denominator = zeros(sz); % [-]
                for j = 1:num_species
                    phi = phi_f(species_list(i).M, species_list(j).M, visc_i{i}, visc_i{j});
                    denominator = denominator + X_i{j} .* phi; 
                end
                k = k + X_i{i} .* k_i{i} ./ denominator;
            end

            % Calculating the effective thermal conductivity [W/(m*K)]
            k = k .* obj.flow.porosity + obj.solid_phase_thermal_conductivity * (1 - obj.flow.porosity);
        end
    
        % Function used to calculate the dynamic viscosity [Pa*s]
        function visc = visc_function(obj)
            species_map = obj.species_manager.species;
            species_list = values(species_map);
            num_species = numel(species_list);

            T = obj.energy.temp;  % Temperature field [K]
            M_mix = obj.species_manager.M_mix; % Mixture molar mass [g/mol]
            Y_tot = obj.species_manager.Y_total; % Sum of mass fractions [-]
            sz = size(T);

            % Preallocating the arrays
            visc = zeros(sz); % Mixture dynamic viscosity [Pa*s]

            X_i = cell(num_species,1); % Molar fractions [-]
            visc_i = cell(num_species, 1); % Pure substances dynamic viscosities [Pa*s]

            for i = 1:num_species
                sp = species_list(i);
                
                % Molar fraction: Y * M_mix / M_i [-]
                X_i{i} = (sp.Y ./ Y_tot) .* M_mix ./ sp.M;
                
                % Pure substance dynamic viscosity [Pa*s]
                visc_i{i} = polyval(sp.visc, T);
            end

            for i = 1:num_species
                denominator = zeros(sz);
                for j = 1:num_species
                    denominator = denominator + X_i{j} * obj.phi_visc(i,j); % [-]
                end
                visc = visc + X_i{i} .* visc_i{i} ./ denominator; % [Pa*s]
            end
        end
    end
end