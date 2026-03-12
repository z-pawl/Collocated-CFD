classdef FlowComponent < IComponent
    properties
        noi (1,1)           % Number of iterations

        % Dependencies
        grid Grid2D {mustBeScalarOrEmpty}


        % Collocated grid approach is used to represent velocity and pressure fields 
        vx (:,:) double             % X - component of velocity [m/s]
        vr (:,:) double             % R - component of velocity [m/s]
        p (:,:) double              % Pressure [Pa]

        % Rhie-chow interpolated velocities [m/s]
        vx_faces (:,:) double
        vr_faces (:,:) double


        % Physical properties
        porosity double             % Porosity [-]
        rho (:,:) double            % Density [kg/m^3]
        visc (:,:) double           % Dynamic viscosity [Pa*s]

        % Source terms
        srcx (:,:) double           % Constant momentum source term in the x direction [kg/(m^2*s^2)]
        src_linx (:,:) double       % Linear momentum source term in the x direction [kg/(m^3*s)]
        srcr (:,:) double           % Constant momentum source term in the r direction [kg/(m^2*s^2)]
        src_linr (:,:) double       % Linear momentum source term in the r direction [kg/(m^3*s)]


        % Functions used to calculate physical properties
        rho_function (1,1) function_handle = @(x) ones(x.grid.sz);
        visc_function (1,1) function_handle = @(x) ones(x.grid.sz);
        
        % Functions used to calculate source terms
        srcx_function (1,1) function_handle = @(x) zeros(x.grid.sz);
        src_linx_function (1,1) function_handle = @(x) zeros(x.grid.sz);
        srcr_function (1,1) function_handle = @(x) zeros(x.grid.sz);
        src_linr_function (1,1) function_handle = @(x) zeros(x.grid.sz);


        % Boundary conditions
        vx_bds Boundaries {mustBeScalarOrEmpty}         % X - component velocity boundary conditions
        vr_bds Boundaries {mustBeScalarOrEmpty}         % R - component velocity boundary conditions
        p_bds Boundaries {mustBeScalarOrEmpty}          % Pressure boundary conditions
        p_corr_bds Boundaries {mustBeScalarOrEmpty}     % Pressure correction boundary conditions

        % Dependencies of the geometric term in the discretized momentum in
        % the radial direction on vr_s, vr_P, vr_n [m]
        % geom_term = k_s * vr_s + k_P * vr_P + k_n * vr_n
        geom_term_dep_s (:,:) double
        geom_term_dep_P (:,:) double
        geom_term_dep_n (:,:) double

        % Solver settings
        tol_v (1,1) double {mustBeNonnegative}          % Convergence criteria for velocity
        tol_cont (1,1) double {mustBeNonnegative}       % Convergence criteria for continuity
        relaxation_factor_v (1,1) double                % Velocity relaxation factor
        relaxation_factor_p (1,1) double                % Pressure relaxation factor
        implicit_relaxation_factor_v (1,1) double = 0.9 % Implicit relaxation factor for velocity to ensure diagonal dominance
        inner_iters (1,1)                               % Number of inner iterations per outer iteration
        solver_iters_v (1,1)                            % Number of velocity solver iterations per inner iteration
        solver_iters_p (1,1)                            % Number of pressure solver iterations per inner iteration
        residual_history_vx (:,1) double                % X component of velocity residual history
        residual_history_vr (:,1) double                % R component of velocity residual history
        residual_history_continuity (:,1) double        % Continuity residual history
        m (1,1) double = 3                              % Number of initial iterations to monitor continuity residual
        max_continuity_residual (1,1) double            % Maximum unscaled continuity residual during first m iterations
    end
    methods
        function obj = FlowComponent(grid, vx, vr, p, porosity, rho_function, visc_function, srcx_function, src_linx_function, srcr_function, src_linr_function, tol_v, tol_cont, relaxation_factor_v, relaxation_factor_p, inner_iters, solver_iters_v, solver_iters_p)
            arguments
                grid (1,1) Grid2D
                vx (:,:) double     % Initial x-component velocity [m/s]
                vr (:,:) double     % Initial r-component velocity field [m/s]
                p  (:,:) double     % Initial pressure field [Pa]
                porosity (:,:) double     % Porosity [-]
                rho_function (1,1) function_handle
                visc_function (1,1) function_handle
                srcx_function (1,1) function_handle
                src_linx_function (1,1) function_handle
                srcr_function (1,1) function_handle
                src_linr_function (1,1) function_handle
                tol_v (1,1) double
                tol_cont (1,1) double
                relaxation_factor_v (1,1) double 
                relaxation_factor_p (1,1) double
                inner_iters (1,1)
                solver_iters_v (1,1)
                solver_iters_p (1,1)
            end

            if ~isequal(size(vx), grid.sz) || ~isequal(size(vr), grid.sz) || ~isequal(size(p), grid.sz)
                error("The initial fields must have the same size as the grid");
            end
            
            % Assigning the properties
            obj.noi = 0;

            obj.grid = grid;

            obj.vx = vx;
            obj.vr = vr;
            obj.p = p;

            obj.vx_faces = zeros(grid.sz + [1 0]);
            obj.vr_faces = zeros(grid.sz + [0 1]);

            obj.porosity = porosity;
            obj.rho_function = rho_function;
            obj.visc_function = visc_function;

            obj.srcx_function = srcx_function;
            obj.src_linx_function = src_linx_function;
            obj.srcr_function = srcr_function;
            obj.src_linr_function = src_linr_function;

            obj.vx_bds = Boundaries;
            obj.vr_bds = Boundaries;
            obj.p_bds = Boundaries;
            obj.p_corr_bds = Boundaries;

            obj.tol_v = tol_v;
            obj.tol_cont = tol_cont;
            obj.relaxation_factor_v = relaxation_factor_v;
            obj.relaxation_factor_p = relaxation_factor_p;
            obj.inner_iters = inner_iters;
            obj.solver_iters_v = solver_iters_v;
            obj.solver_iters_p = solver_iters_p;
            obj.residual_history_vx = zeros(10000,1);
            obj.residual_history_vr = zeros(10000,1);
            obj.residual_history_continuity = zeros(10000,1);

            % Geometric term dependence
            % The radial velocity profile is approximated using a quadratic
            % function fitted using the radial velocity in the cell center
            % and north and south faces
            % vr(r) = a * r^2 + b * r + c
            % vr(xs) = vr(rP - dr /2) = vr_s; vr(rP) = vr_P; vr(rn) = vr(rP + dr /2) = vr_n;
            % Then the geometric term vr/r^2 is integrated over a control
            % volume (phi=0..2pi, r=rs..rn=rP-dr/2..rP+dr/2, x=xw..xe=xP-dx/2..xP+dx/2

            % Following expressions are then obtained
            k = 2 * pi * obj.grid.dx ./ (obj.grid.dr .^ 2); % [1/m]
            log_term = log((obj.grid.cent_pos_r + obj.grid.dr / 2) / (obj.grid.cent_pos_r - obj.grid.dr / 2)); % [-]

            obj.geom_term_dep_s = k .* (log_term .* (2 * obj.grid.cent_pos_r .^ 2 + obj.grid.dr .* obj.grid.cent_pos_r) ...
                - obj.grid.dr .^ 2 - 2 * obj.grid.dr .* obj.grid.cent_pos_r); % [m]

            obj.geom_term_dep_P = k .* (log_term .* (obj.grid.dr .^ 2 - 4 * obj.grid.cent_pos_r .^ 2) ...
                + 4 * obj.grid.dr .* obj.grid.cent_pos_r); % [m]

            obj.geom_term_dep_n = k .* (log_term .* (2 * obj.grid.cent_pos_r .^ 2 - obj.grid.dr .* obj.grid.cent_pos_r) ...
                + obj.grid.dr .^ 2 - 2 * obj.grid.dr .* obj.grid.cent_pos_r); % [m]
        end

        % Updates fluid properties and the source terms
        function update_properties(obj)
            obj.rho = obj.rho_function();
            obj.visc = obj.visc_function();

            obj.srcx = obj.srcx_function();
            obj.src_linx = obj.src_linx_function();
            obj.srcr = obj.srcr_function();
            obj.src_linr = obj.src_linr_function();
        end

        % Updates face velocities, face velocities are calculated using Rhie-Chow interpolation
        function update_face_velocities(obj, vx, vr, aP_vx, aP_vr)
            [rhie_chow_vx, rhie_chow_vr] = rhie_chow_interpolation(obj.grid, vx, vr, obj.p, obj.vx_bds, obj.vr_bds, obj.p_bds, aP_vx, aP_vr);
            
            obj.vx_faces = rhie_chow_vx; % [m/s]
            obj.vr_faces = rhie_chow_vr; % [m/s]
        end

        % Converts pressure boundary conditions to pressure correction boundary conditions
        function convert_pressure_bd_conditions(obj)
            obj.p_corr_bds = Boundaries;

            for i = 1:size(obj.p_bds.boundaries,1)
                obj.p_corr_bds.add_boundary(obj.p_bds.boundaries(i).convert_into_correction_boundary())
            end
        end

        function [coeff_vx, coeff_vr] = get_coefficients_v(obj)
            % Coefficients used to calculate the physical properties at faces
            [lin_int_coeffs_x, lin_int_coeffs_r] = linear_interpolation_scheme(obj.grid);
            [lin_int_coeffs_x, lin_int_coeffs_r] = obj.grid.domain_boundary.apply_boundary_condition_value(lin_int_coeffs_x, lin_int_coeffs_r);
            
            % Properties interpolated to faces
            % Density [kg/m^3]
            [rho_x, rho_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.rho);
            % Dynamic viscosity [Pa*s]
            [visc_x, visc_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.visc);
            % Porosity [-]
            [porosity_x, porosity_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.porosity);

            clear lin_int_coeffs_x lin_int_coeffs_r

            % F = rho * v * A / porosity ^ 2 [kg/s]
            Fx = rho_x .* obj.vx_faces .* obj.grid.face_area_x ./ porosity_x .^ 2;
            Fr = rho_r .* obj.vr_faces .* obj.grid.face_area_r ./ porosity_r .^ 2;

            % D = visc * A / porosity [kg*m/s]
            Dx = visc_x .* obj.grid.face_area_x ./ porosity_x;
            Dr = visc_r .* obj.grid.face_area_r ./ porosity_r;

            clear rho_x rho_r visc_x visc_r;


            % Calculating velocities coefficients
            coeff_vx = obj.get_coefficients_vx(Fx, Fr, Dx, Dr);
            coeff_vr = obj.get_coefficients_vr(Fx, Fr, Dx, Dr);
            clear Fx Fr Dx Dr;


            % Calculating the volume integral of pressure gradient [N]=[kg*m/s^2]

            % Coefficients for pressure interpolation
            [coeff_px, coeff_pr] = linear_interpolation_scheme(obj.grid);
            [coeff_px, coeff_pr] = obj.p_bds.apply_boundary_condition_value(coeff_px, coeff_pr);

            % Pressure values at faces [Pa]
            [px, pr] = evaluate_faces(coeff_px, coeff_pr, obj.p);

            % Adding the volume integral of pressure gradient as source terms
            coeff_vx(:,:,6) = coeff_vx(:,:,6) + obj.grid.volume ./ obj.grid.dx .* (px(1:end-1,:) - px(2:end,:));
            coeff_vr(:,:,6) = coeff_vr(:,:,6) + obj.grid.volume ./ obj.grid.dr .* (pr(:,1:end-1) - pr(:,2:end));

            % Adding implicit relaxation to ensure the diagonal dominance
            coeff_vx = implicit_relaxation(coeff_vx, obj.vx, obj.implicit_relaxation_factor_v);
            coeff_vr = implicit_relaxation(coeff_vr, obj.vr, obj.implicit_relaxation_factor_v);
        end

        function coeff_vx = get_coefficients_vx(obj, Fx, Fr, Dx, Dr)
            % Values and normal derivatives of velocity at faces

            % Coefficients for normal derivatives
            [coeffs_x_n_der_vx, coeffs_r_n_der_vx] = central_differencing_scheme(obj.grid);
            [coeffs_x_n_der_vx, coeffs_r_n_der_vx] = obj.vx_bds.apply_boundary_condition_normal_derivative(coeffs_x_n_der_vx, coeffs_r_n_der_vx);

            % Values of normal derivatives [1/s]
            [n_der_vx_x, n_der_vx_r] = evaluate_faces(coeffs_x_n_der_vx, coeffs_r_n_der_vx, obj.vx);

            % Coefficients for the upwind scheme
            [coeffs_x_upwind, coeffs_r_upwind] = upwind_scheme(obj.grid, obj.vx_faces, obj.vr_faces);

            % Coefficients for the TVD scheme
            [coeffs_x_TVD, coeffs_r_TVD] = TVD_scheme(obj.grid, obj.vx_faces, obj.vr_faces, @van_leer_flux_limiter, n_der_vx_x, n_der_vx_r);

            % Coefficients for the deferred scheme
            [coeffs_x_vx_deferred, coeffs_r_vx_deferred] = deferred_correction_face(coeffs_x_upwind, coeffs_r_upwind, coeffs_x_TVD, coeffs_r_TVD, obj.vx);
            [coeffs_x_vx_deferred, coeffs_r_vx_deferred] = obj.vx_bds.apply_boundary_condition_value(coeffs_x_vx_deferred, coeffs_r_vx_deferred);

            clear n_der_vx_x n_der_vx_r coeffs_x_upwind coeffs_r_upwind coeffs_x_TVD coeffs_r_TVD;

            coeff_vx = assemble_coeff_array(Fx, Fr, Dx, Dr, coeffs_x_vx_deferred, coeffs_r_vx_deferred, coeffs_x_n_der_vx, coeffs_r_n_der_vx, obj.srcx, obj.src_linx, obj.grid.volume);
        end

        function coeff_vr = get_coefficients_vr(obj, Fx, Fr, Dx, Dr)
            % Values and normal derivatives of velocity at faces

            % Coefficients for normal derivatives
            [coeffs_x_n_der_vr, coeffs_r_n_der_vr] = central_differencing_scheme(obj.grid);
            [coeffs_x_n_der_vr, coeffs_r_n_der_vr] = obj.vr_bds.apply_boundary_condition_normal_derivative(coeffs_x_n_der_vr, coeffs_r_n_der_vr);

            % Values of normal derivatives [1/s]
            [n_der_vr_x, n_der_vr_r] = evaluate_faces(coeffs_x_n_der_vr, coeffs_r_n_der_vr, obj.vr);

            % Coefficients for the upwind scheme
            [coeffs_x_upwind, coeffs_r_upwind] = upwind_scheme(obj.grid, obj.vx_faces, obj.vr_faces);

            % Coefficients for the TVD scheme
            [coeffs_x_TVD, coeffs_r_TVD] = TVD_scheme(obj.grid, obj.vx_faces, obj.vr_faces, @van_leer_flux_limiter, n_der_vr_x, n_der_vr_r);

            % Coefficients for the deferred scheme
            [coeffs_x_vr_deferred, coeffs_r_vr_deferred] = deferred_correction_face(coeffs_x_upwind, coeffs_r_upwind, coeffs_x_TVD, coeffs_r_TVD, obj.vr);
            [coeffs_x_vr_deferred, coeffs_r_vr_deferred] = obj.vr_bds.apply_boundary_condition_value(coeffs_x_vr_deferred, coeffs_r_vr_deferred);

            clear n_der_vr_x n_der_vr_r coeffs_x_upwind coeffs_r_upwind coeffs_x_TVD coeffs_r_TVD;

        
            coeff_vr = assemble_coeff_array(Fx, Fr, Dx, Dr, coeffs_x_vr_deferred, coeffs_r_vr_deferred, coeffs_x_n_der_vr, coeffs_r_n_der_vr, obj.srcr, obj.src_linr, obj.grid.volume);
        
            % Additional geometric term
            [~, vr_at_r_faces] = evaluate_faces(coeffs_x_vr_deferred, coeffs_r_vr_deferred, obj.vr); % [m/s]
            geom_term = obj.visc ./ obj.porosity ...
                .* (obj.geom_term_dep_s .* vr_at_r_faces(:,1:end-1) ...
                + obj.geom_term_dep_n .* vr_at_r_faces(:,2:end) ...
                + obj.geom_term_dep_P .* obj.vr); % [N]=[kg*m/s^2]

            coeff_vr(:,:,6) = coeff_vr(:,:,6) - geom_term;
        end

        function coeff_p_corr = get_coefficients_p_corr(obj, coeff_vx_aP, coeff_vr_aP)
            % Coefficients used to calculate the linear interpolation
            [lin_int_coeffs_x, lin_int_coeffs_r] = linear_interpolation_scheme(obj.grid);
            [lin_int_coeffs_x, lin_int_coeffs_r] = obj.grid.domain_boundary.apply_boundary_condition_value(lin_int_coeffs_x, lin_int_coeffs_r);

            % Interpolated quantities
            % Density [kg/m^3]
            [rho_x, rho_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.rho);
            % Df [m^3*s/kg]
            [Df_x, ~] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.grid.volume ./ coeff_vx_aP);
            [~, Df_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.grid.volume ./ coeff_vr_aP);

            clear lin_int_coeffs_x lin_int_coeffs_x coeff_vx_aP coeff_vr_aP

            % Calculating the coefficients used in the discretized equations
            c_x = rho_x .* obj.grid.face_area_x .* Df_x; % [m^2*s]
            c_r = rho_r .* obj.grid.face_area_r .* Df_r; % [m^2*s]
            Fx = rho_x .* obj.vx_faces .* obj.grid.face_area_x; % [kg/s]
            Fr = rho_r .* obj.vr_faces .* obj.grid.face_area_r; % [kg/s]
            clear rho_x rho_r coeff_vx_aP_f coeff_vr_aP_f Df_x Df_r

            % Calculating coefficients for normal derivatives at faces for pressure
            [coeff_der_p_x, coeff_der_p_r] = central_differencing_scheme(obj.grid);
            [coeff_der_p_x, coeff_der_p_r] = obj.p_corr_bds.apply_boundary_condition_normal_derivative(coeff_der_p_x, coeff_der_p_r);


            % Preallocating the array
            coeff_p_corr = zeros([obj.grid.sz 6]);

            % aP [m*s]
            coeff_p_corr(:,:,1) = -c_x(1:end-1,:) .* coeff_der_p_x(1:end-1,:,2) + c_x(2:end,:) .* coeff_der_p_x(2:end,:,1) ...
                + -c_r(:,1:end-1) .* coeff_der_p_r(:,1:end-1,2) + c_r(:,2:end) .* coeff_der_p_r(:,2:end,1);
            % aW [m*s]
            coeff_p_corr(:,:,2) = c_x(1:end-1,:) .* coeff_der_p_x(1:end-1,:,1);
            % aE [m*s]
            coeff_p_corr(:,:,3) = -c_x(2:end,:) .* coeff_der_p_x(2:end,:,2);
            % aS [m*s]
            coeff_p_corr(:,:,4) = c_r(:,1:end-1) .* coeff_der_p_r(:,1:end-1,1);
            % aN [m*s]
            coeff_p_corr(:,:,5) = -c_r(:,2:end) .* coeff_der_p_r(:,2:end,2);
            % b [kg/s]
            coeff_p_corr(:,:,6) = Fx(2:end,:) - Fx(1:end-1,:) + Fr(:,2:end) - Fr(:,1:end-1) ...
                + c_x(1:end-1,:) .* coeff_der_p_x(1:end-1,:,3) - c_x(2:end,:) .* coeff_der_p_x(2:end,:,3) ...
                + c_r(:,1:end-1) .* coeff_der_p_r(:,1:end-1,3) - c_r(:,2:end) .* coeff_der_p_r(:,2:end,3);
        end

        function converged = iterate(obj)
            obj.noi = obj.noi + 1;

            % Inner iterations
            for i = 1:obj.inner_iters
                % Updating properties and calculating coefficients of velocity discretized equations
                obj.update_properties();
                [coeff_vx, coeff_vr] = obj.get_coefficients_v();

                % Calculating the intermediate velocities [m/s]
                vx_star=solve(coeff_vx, obj.vx, obj.solver_iters_v, randi([1 4], [4 1]), obj.relaxation_factor_v);
                vr_star=solve(coeff_vr, obj.vr, obj.solver_iters_v, randi([1 4], [4 1]), obj.relaxation_factor_v);

                % Calculating the unrelaxed central coefficients used to
                % calculate rhie chow correction and the pressure correction [kg/s]
                coeff_vx_unrelaxed = coeff_vx(:,:,1) * obj.implicit_relaxation_factor_v;
                coeff_vr_unrelaxed = coeff_vr(:,:,1) * obj.implicit_relaxation_factor_v;

                % Calculating the face velocities using the intermediate velocities
                obj.update_face_velocities(vx_star, vr_star, coeff_vx_unrelaxed, coeff_vr_unrelaxed);

                % Calculating the coefficients of the pressure correction equation and solving it
                coeff_p_corr = obj.get_coefficients_p_corr(coeff_vx_unrelaxed, coeff_vr_unrelaxed);

                p_corr = solve(coeff_p_corr, zeros(obj.grid.sz), obj.solver_iters_p, [1;3;2;4], 1); % [Pa]

                % Calculating the pressure force from the pressure correction and using it to correct the velocities
                % Coefficients for pressure interpolation
                [coeff_px, coeff_pr] = linear_interpolation_scheme(obj.grid);
                [coeff_px, coeff_pr] = obj.p_corr_bds.apply_boundary_condition_value(coeff_px, coeff_pr);
                % Pressure correction values at faces [Pa]
                [p_corr_x, p_corr_r] = evaluate_faces(coeff_px, coeff_pr, p_corr);
                % Adding the volume integral of pressure gradient as source terms
                p_force_x = obj.grid.volume ./ obj.grid.dx .* (p_corr_x(1:end-1,:) - p_corr_x(2:end,:)); % [N]=[kg*m/s^2]
                p_force_r = obj.grid.volume ./ obj.grid.dr .* (p_corr_r(:,1:end-1) - p_corr_r(:,2:end)); % [N]=[kg*m/s^2]

                vx_star = vx_star + p_force_x ./ coeff_vx(:,:,1) * obj.relaxation_factor_v; % [m/s]
                vr_star = vr_star + p_force_r ./ coeff_vr(:,:,1) * obj.relaxation_factor_v; % [m/s]

                % Updating the fields using the relaxation factor
                obj.p = obj.p + p_corr * obj.relaxation_factor_p; % [Pa]
                obj.vx = vx_star; % [m/s]
                obj.vr = vr_star; % [m/s]

                % Updating the face velocities and calculating the face velocities using new velocities
                obj.update_properties();
                [coeff_vx, coeff_vr] = obj.get_coefficients_v();
                coeff_vx_unrelaxed = coeff_vx(:,:,1) * obj.implicit_relaxation_factor_v; % [kg/s]
                coeff_vr_unrelaxed = coeff_vr(:,:,1) * obj.implicit_relaxation_factor_v; % [kg/s]
                obj.update_face_velocities(obj.vx, obj.vr, coeff_vx_unrelaxed, coeff_vr_unrelaxed);
            end

            % Updating properties
            obj.update_properties();

            % Calculating the residual and checking convergence
            [coeff_vx, coeff_vr] = obj.get_coefficients_v();
            coeff_vx_unrelaxed = coeff_vx(:,:,1) * obj.implicit_relaxation_factor_v; % [kg/s]
            coeff_vr_unrelaxed = coeff_vr(:,:,1) * obj.implicit_relaxation_factor_v; % [kg/s]
            obj.update_face_velocities(obj.vx, obj.vr, coeff_vx_unrelaxed, coeff_vr_unrelaxed);

            res_vx = calculate_residual(coeff_vx, obj.vx);
            res_vr = calculate_residual(coeff_vr, obj.vr);

            % The density has to interpolated in order to calculate the mass
            % creation rate in each cell
            [coeff_lin_x, coeff_lin_r] = linear_interpolation_scheme(obj.grid);
            [coeff_lin_x, coeff_lin_r] = obj.grid.domain_boundary.apply_boundary_condition_value(coeff_lin_x, coeff_lin_r);
            [rho_x, rho_r] = evaluate_faces(coeff_lin_x, coeff_lin_r, obj.rho);
            mass_flux_x = rho_x .* obj.grid.face_area_x .* obj.vx_faces;
            mass_flux_r = rho_r .* obj.grid.face_area_r .* obj.vr_faces;

            res_continuity = norm((mass_flux_x(2:end,:) - mass_flux_x(1:end-1,:)) + (mass_flux_r(:,2:end) - mass_flux_r(:,1:end-1)), 1);
            if obj.noi <= obj.m
                obj.max_continuity_residual = max([obj.max_continuity_residual res_continuity]);
            end
            res_continuity = res_continuity / obj.max_continuity_residual;

            obj.residual_history_vx(obj.noi) = res_vx;
            obj.residual_history_vr(obj.noi) = res_vr;
            obj.residual_history_continuity(obj.noi) = res_continuity;
            
            converged = (res_vx < obj.tol_v) && (res_vr < obj.tol_v) && (res_continuity < obj.tol_cont);
        end
    end
end