classdef SpeciesComponent < handle
    properties
        name (1,1) string           % Species name

        % Dependencies
        grid Grid2D {mustBeScalarOrEmpty}
        flow FlowComponent {mustBeScalarOrEmpty}
        species_manager SpeciesManagerComponent {mustBeScalarOrEmpty}

        % Physical properties
        Y (:,:) double              % Mass fraction [-]
        D (:,:) double              % Effective mass diffusivity in mixture [m^2/s]
        p_partial (:,:) double      % Partial pressure [Pa]

        % Pure substance properties
        Cp double                   % Molar heat capacity [J/(mol*K)] - polynomial of temperature [K]
        visc double                 % Dynamic viscosity [Pa*s] - polynomial of temperature [K]
        k double                    % Thermal conductivity [W/(m*K)] - polynomial of temperature [K]
        M double                    % Molar mass [g/mol]
        V_D double                  % Diffusion volume [?]

        % Source terms
        src (:,:) double            % Constant source term [kg/(m^3*s)]
        src_lin (:,:) double        % Linear source term [kg/(m^3*s)]

        % Boundary conditions
        species_bds Boundaries {mustBeScalarOrEmpty}    % Species boundary conditions

        % Solver settings
        tol (1,1) double {mustBeNonnegative}            % Convergence criteria [-]
        relaxation_factor (1,1) double                  % Relaxation factor
        solver_iters (1,1)                              % Number of solver iterations per inner iteration

        aP (:,:) double                                 % Central coefficient of the discretized equation [kg/s]
        Jx (:,:) double                                 % Diffusive flux in the x direction [kg/s]
        Jr (:,:) double                                 % Diffusive flux in the r direction [kg/s]
    end
    methods
        function obj = SpeciesComponent(name, grid, flow, species_manager, Y, Cp, visc, k, M, V_D, tol, relaxation_factor, solver_iters)
            arguments
                name (1,1) string               % Name of the species
                grid (1,1) Grid2D
                flow (1,1) FlowComponent
                species_manager (1,1) SpeciesManagerComponent
                Y (:,:) double                  % Initial mass fraction field [-]
                Cp (1,:) double                 % Polynomial describing the heat capacity of the pure substance [J/(mol*K)] as a function of temperature [K]
                visc (1,:) double               % Polynomial describing the dynamic viscosity of the pure substance [Pa*s] as a function of temperature [K]
                k (1,:) double                  % % Polynomial describing the heat conductivity of the pure substance [W/(m*K)] as a function of temperature [K]
                M (1,1) double                  % Molar mass of the pure substance [g/mol]
                V_D (1,1) double                % Diffusion volume [?]
                tol (1,1) double                % Convergence criteria [-]
                relaxation_factor (1,1) double  % Relaxation factor
                solver_iters (1,1)              % Number of solver iterations per inner iteration
            end

            if ~isequal(size(Y), grid.sz)
                error("The initial fields must have the same size as the grid");
            end
            
            % Assigning the properties
            obj.name = name;

            obj.grid = grid;
            obj.flow = flow;
            obj.species_manager = species_manager;

            obj.Y = Y;
            obj.D = ones(grid.sz);
            obj.p_partial = zeros(grid.sz);

            obj.Cp = Cp;
            obj.visc = visc;
            obj.k = k;
            obj.M = M;
            obj.V_D = V_D;

            obj.src = zeros(grid.sz);
            obj.src_lin = zeros(grid.sz);

            obj.species_bds = Boundaries;

            obj.tol = tol;
            obj.relaxation_factor = relaxation_factor;
            obj.solver_iters = solver_iters;

            obj.aP = ones(grid.sz);
            obj.Jx = zeros(grid.sz + [1 0]);
            obj.Jr = zeros(grid.sz + [0 1]);

            species_manager.add_species(name, obj);
        end

        function coeff = get_coefficients(obj)
            % Coefficients used to calculate the physical properties at faces
            [lin_int_coeffs_x, lin_int_coeffs_r] = linear_interpolation_scheme(obj.grid);
            [lin_int_coeffs_x, lin_int_coeffs_r] = obj.grid.domain_boundary.apply_boundary_condition_value(lin_int_coeffs_x, lin_int_coeffs_r);
            
            % Properties interpolated to faces
            % Density [kg/m^3]
            [rho_x, rho_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.flow.rho);
            % Diffusivity [m^2/s]
            [diff_x, diff_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.D);

            % Values and normal derivatives of temperature at faces
            % Coefficients for normal derivatives
            [coeffs_x_n_der, coeffs_r_n_der] = central_differencing_scheme(obj.grid);
            [coeffs_x_n_der, coeffs_r_n_der] = obj.species_bds.apply_boundary_condition_normal_derivative(coeffs_x_n_der, coeffs_r_n_der);

            % Values of normal derivatives [1/m]
            [n_der_x, n_der_r] = evaluate_faces(coeffs_x_n_der, coeffs_r_n_der, obj.Y);

            % Coefficients for the upwind scheme
            [coeffs_x_upwind, coeffs_r_upwind] = upwind_scheme(obj.grid, obj.flow.vx_faces, obj.flow.vr_faces);
            
            % Coefficients for the TVD scheme
            [coeffs_x_TVD, coeffs_r_TVD] = TVD_scheme(obj.grid, obj.flow.vx_faces, obj.flow.vr_faces, @van_leer_flux_limiter, n_der_x, n_der_r);

            % Coefficients for the deferred scheme
            [coeffs_x_deferred, coeffs_r_deferred] = deferred_correction_face(coeffs_x_upwind, coeffs_r_upwind, coeffs_x_TVD, coeffs_r_TVD, obj.Y);
            [coeffs_x_deferred, coeffs_r_deferred] = obj.species_bds.apply_boundary_condition_value(coeffs_x_deferred, coeffs_r_deferred);

            % F = rho * v * A [kg/s]
            Fx = rho_x .* obj.flow.vx_faces .* obj.grid.face_area_x;
            Fr = rho_r .* obj.flow.vr_faces .* obj.grid.face_area_r;

            % D = rho * diff * A [kg*m/s]
            Dx = rho_x .* diff_x .* obj.grid.face_area_x;
            Dr = rho_r .* diff_r .* obj.grid.face_area_r;

            % The sum of diffusive fluxes for all species is subtracted
            % from the mass flux in order to apply mass correction
            % J* = J - Y_i * J_sum; J - diffusive flux [kg/s]
            coeff = assemble_coeff_array(Fx - obj.species_manager.Jx_sum, Fr - obj.species_manager.Jr_sum, Dx, Dr, coeffs_x_deferred, coeffs_r_deferred, coeffs_x_n_der, coeffs_r_n_der, obj.src, obj.src_lin, obj.grid.volume);
        end

        function inner_iteration(obj)
            % This is a single inner iteration as the outer iteration is
            % handled by the SpeciesManagerComponent class

            % Calculating coefficients of discretized equations
            coeff = obj.get_coefficients();

            % Solving the system of equations and updating the field using the relaxation factor
            obj.Y = obj.Y + obj.relaxation_factor * (solve_mex(coeff, obj.Y, obj.solver_iters, [1; 3; 2; 4], 1) - obj.Y);
            

            % Updating the aP, Jx and Jr values used by the species manager

            obj.aP = coeff(:,:,1);

            % Coefficients used to calculate the physical properties at faces
            [lin_int_coeffs_x, lin_int_coeffs_r] = linear_interpolation_scheme(obj.grid);
            [lin_int_coeffs_x, lin_int_coeffs_r] = obj.grid.domain_boundary.apply_boundary_condition_value(lin_int_coeffs_x, lin_int_coeffs_r);
            
            % Properties interpolated to faces
            % Density [kg/m^3]
            [rho_x, rho_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.flow.rho);
            % Diffusivity [m^2/s]
            [diff_x, diff_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.D);

            % Values and normal derivatives of temperature at faces
            % Coefficients for normal derivatives
            [coeffs_x_n_der, coeffs_r_n_der] = central_differencing_scheme(obj.grid);
            [coeffs_x_n_der, coeffs_r_n_der] = obj.species_bds.apply_boundary_condition_normal_derivative(coeffs_x_n_der, coeffs_r_n_der);

            % Values of normal derivatives [1/m]
            [n_der_x, n_der_r] = evaluate_faces(coeffs_x_n_der, coeffs_r_n_der, obj.Y);

            % D = rho * diff * A [kg*m/s]
            Dx = rho_x .* diff_x .* obj.grid.face_area_x;
            Dr = rho_r .* diff_r .* obj.grid.face_area_r;

            % Diffusive fluxes along faces [kg/s]
            obj.Jx = -Dx .* n_der_x;
            obj.Jr = -Dr .* n_der_r;
        end

        function res = calculate_scaled_residual(obj)
            res = calculate_residual(obj.get_coefficients(), obj.Y);
        end
    end
end