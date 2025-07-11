classdef ThermophysicalProperties < handle
    properties
        % Dependencies
        grid Grid2D {mustBeScalarOrEmpty}
        flow FlowComponent {mustBeScalarOrEmpty}
        energy EnergyComponent {mustBeScalarOrEmpty}
        species_manager SpeciesManagerComponent {mustBeScalarOrEmpty}

        D_AB double             % Binary diffusion coefficients (without the part depending on temperature and pressure)
        phi_visc double         % Coefficients used to calculate the mixture viscosity
    end
    methods
        function obj = ThermophysicalProperties(grid, flow, energy, species_manager)
            arguments
                grid (1,1) Grid2D
                flow (1,1) FlowComponent
                energy (1,1) EnergyComponent
                species_manager (1,1) SpeciesManagerComponent
            end

            obj.grid = grid;
            obj.flow = flow;
            obj.energy = energy;
            obj.species_manager = species_manager;

            obj.update_values();
        end
        
        function update_values(obj)
            obj.calculate_D_AB();
            obj.calculate_phi_visc();
        end

        function calculate_D_AB(obj)
            C = 0.00143;  % Empirical constant
            species_list = obj.species_manager.species;
            species_list = values(species_list);
            num_species = numel(species_list);
        
            % Preallocate diffusivity matrix
            D_AB = zeros(num_species);
        
            % Precompute properties
            M = zeros(num_species, 1);      % Molar masses
            V_D = zeros(num_species, 1);    % Diffusion volumes
        
            for i = 1:num_species
                M(i) = species_list(i).M;
                V_D(i) = species_list(i).V_D;
            end
        
            % Compute binary diffusion coefficients
            for i = 1:num_species
                for j = 1:num_species
                    if i == j
                        D_AB(i, j) = Inf;
                    else
                        M_AB = 2 / (1 / M(i) + 1 / M(j));
                        sigma = (V_D(i)^(1/3) + V_D(j)^(1/3))^2;
                        D_AB(i, j) = C / (sqrt(M_AB) * sigma);
                    end
                end
            end
        
            obj.D_AB = D_AB;
        end

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
            M_mix = obj.species_manager.M_mix;
            molar_fractions = cell(num_species,1);
            
            for i = 1:num_species
                sp = species_list(i);
                molar_fractions{i} = sp.Y .* M_mix ./ sp.M;
            end
        
            % Allocate output map for diffusivities
            D = containers.Map();
        
            % Common parameters
            T = obj.energy.temp;     % Temperature field
            P = obj.flow.p;          % Pressure field
        
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
                D(name_i) = D_i;
            end
        end

        function rho = rho_function(obj)
            % The density is calculated using the ideal gas law
            rho = obj.flow.p .* obj.species_manager.M_mix ./ (obj.energy.temp * obj.species_manager.R);
            

            % g/m^3 to kg/m^3
            unit_conversion = 1/1000;
            rho = rho * unit_conversion;
        end
    
        function cp = cp_function(obj)
            species_map = obj.species_manager.species;
            species_list = values(species_map);
            num_species = numel(species_list);
            
            T = obj.energy.temp;  % Temperature field
            M_mix = obj.species_manager.M_mix;
        
            % Preallocate arrays
            cp = zeros(size(T));
            
            for i = 1:num_species
                sp = species_list(i);
                
                % Molar fraction: Y * M_mix / M_i
                X_i = sp.Y .* M_mix ./ sp.M;
                
                % Heat capacity
                cp_i = polyval(sp.Cp, T);
                
                cp = cp + X_i .* cp_i;
            end
        
            % Convert from molar to specific heat capacity
            JGK_TO_JKGK = 1 / 1000;
            cp = cp .* M_mix * JGK_TO_JKGK;
        end
    
        function k = k_function(obj)
            % Nested helper function for phi
            function phi = phi_f(Mi, Mj, visc_i, visc_j) 
                M_ratio = Mi / Mj;
                visc_ratio = visc_i ./ visc_j;

                phi = (1 + sqrt(visc_ratio) * M_ratio ^ (-0.25)) .^ 2 / (8 * (1 + M_ratio)) ^ (1/2);
            end

            species_map = obj.species_manager.species;
            species_list = values(species_map);
            num_species = numel(species_list);

            T = obj.energy.temp;  % Temperature field
            M_mix = obj.species_manager.M_mix;
            sz = size(T);

            % Preallocating the arrays
            k = zeros(sz); % Mixture thermal conductivity

            X_i = cell(num_species,1); % Molar fractions
            visc_i = cell(num_species, 1); % Dynamic viscosities
            k_i = cell(num_species, 1); % Thermal conductivities

            for i = 1:num_species
                sp = species_list(i);
                
                % Molar fraction: Y * M_mix / M_i
                X_i{i} = sp.Y .* M_mix ./ sp.M;
                
                % Dynamic viscosity
                visc_i{i} = polyval(sp.visc, T);

                % Thermal conductivity
                k_i{i} = polyval(sp.k, T);
            end

            for i = 1:num_species
                denominator = zeros(sz);
                for j = 1:num_species
                    phi = phi_f(species_list(i).M, species_list(j).M, visc_i{i}, visc_i{j});
                    denominator = denominator + X_i{j} .* phi;
                end
                k = k + X_i{i} .* k_i{i} ./ denominator;
            end
        end
    
        function visc = visc_function(obj)
            species_map = obj.species_manager.species;
            species_list = values(species_map);
            num_species = numel(species_list);

            T = obj.energy.temp;  % Temperature field
            M_mix = obj.species_manager.M_mix;
            sz = size(T);

            % Preallocating the arrays
            visc = zeros(sz); % Mixture dynamic viscosity

            X_i = cell(num_species,1); % Molar fractions
            visc_i = cell(num_species, 1); % Dynamic viscosities

            for i = 1:num_species
                sp = species_list(i);
                
                % Molar fraction: Y * M_mix / M_i
                X_i{i} = sp.Y .* M_mix ./ sp.M;
                
                % Dynamic viscosity
                visc_i{i} = polyval(sp.visc, T);
            end

            for i = 1:num_species
                denominator = zeros(sz);
                for j = 1:num_species
                    denominator = denominator + X_i{j} * obj.phi_visc(i,j);
                end
                visc = visc + X_i{i} .* visc_i{i} ./ denominator;
            end
        end
    end
end