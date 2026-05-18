classdef NonUniformConvectiveBoundary < IBoundaryType
    properties
        grid Grid2D {mustBeScalarOrEmpty}
        energy_component EnergyComponent {mustBeScalarOrEmpty}

        % Linearized indices
        lin_ind_x_p1 (:,1) double % x face, page 1 (i,j,1)
        lin_ind_x_p2 (:,1) double % x face, page 2 (i,j,2)
        lin_ind_x_p3 (:,1) double % x face, page 3 (i,j,3)

        lin_ind_r_p1 (:,1) double % r face, page 1 (i,j,1)
        lin_ind_r_p2 (:,1) double % r face, page 2 (i,j,2)
        lin_ind_r_p3 (:,1) double % r face, page 3 (i,j,3)

        % Coefficient values
        % X faces
        heat_trans_coeff_x (:,1) double % Heat transfer coefficient [W/(m^2*K)]
        environment_temp_x (:,1) double % Environment temperature [K]
        dx_inv_times_2 (:,1) double % Inverse of the distance to the boundary [1/m]
        orientation_x (:,1) double % Orientation - points inwards [-]
        cell_ind_x (:,1) double % Linearized index of the boundary cell

        % R faces
        heat_trans_coeff_r (:,1) double % Heat transfer coefficient [W/(m^2*K)]
        environment_temp_r (:,1) double % Environment temperature [K]
        dr_inv_times_2 (:,1) double % Inverse of the distance to the boundary [1/m]
        orientation_r (:,1) double % Orientation - points inwards [-]
        cell_ind_r (:,1) double % Linearized index of the boundary cell
    end
    methods
        function obj = NonUniformConvectiveBoundary(name, grid, energy_component)
            arguments
                name (1,1) string
                grid (1,1) Grid2D
                energy_component (1,1) EnergyComponent
            end
            obj.name = name;
            obj.grid = grid;
            obj.energy_component = energy_component;
        end

        function add_boundary_face(obj, direction, face_index, orientation, heat_transfer_coefficient, environment_temperature)
            arguments
                obj
                direction (1,1) string
                face_index (1,2) double
                orientation (1,1) double
                heat_transfer_coefficient (1,1) double
                environment_temperature (1,1) double
            end

            if orientation ~= -1 && orientation ~= 1
                error("The orientation has to be either -1 or 1");
            end

            switch direction
                case "x"
                    if (1 > face_index(1) || face_index(1) > obj.grid.sz(1)+1) || (1 > face_index(2) || face_index(2) > obj.grid.sz(2))
                        error("Provided face index is out of bounds");
                    end

                    % Linearizing the index
                    new_ind = numel(obj.lin_ind_x_p1) + 1;
                    sz = obj.grid.sz + [1 0];
                    lin_ind = sub2ind(sz, face_index(1), face_index(2));
                    obj.lin_ind_x_p1(new_ind) = lin_ind;
                    obj.lin_ind_x_p2(new_ind) = lin_ind + prod(sz);
                    obj.lin_ind_x_p3(new_ind) = lin_ind + 2 * prod(sz);

                    % Coefficients
                    obj.heat_trans_coeff_x(new_ind) = heat_transfer_coefficient;
                    obj.environment_temp_x(new_ind) = environment_temperature;
                    obj.dx_inv_times_2(new_ind) = 2 / obj.grid.dx(face_index(1) + (orientation - 1) / 2);
                    obj.orientation_x(new_ind) = orientation;
                    obj.cell_ind_x(new_ind) = sub2ind(obj.grid.sz, face_index(1) + (orientation - 1) / 2, face_index(2));
                case "r"
                    if (1 > face_index(1) || face_index(1) > obj.grid.sz(1)) || (1 > face_index(2) || face_index(2) > obj.grid.sz(2)+1)
                        error("Provided face index is out of bounds");
                    end
                    
                    % Linearizing the index
                    new_ind = numel(obj.lin_ind_r_p1) + 1;
                    sz = obj.grid.sz + [0 1];
                    lin_ind = sub2ind(sz, face_index(1), face_index(2));
                    obj.lin_ind_r_p1(new_ind) = lin_ind;
                    obj.lin_ind_r_p2(new_ind) = lin_ind + prod(sz);
                    obj.lin_ind_r_p3(new_ind) = lin_ind + 2 * prod(sz);

                    % Coefficients
                    obj.heat_trans_coeff_r(new_ind) = heat_transfer_coefficient;
                    obj.environment_temp_r(new_ind) = environment_temperature;
                    obj.dr_inv_times_2(new_ind) = 2 / obj.grid.dr(face_index(2) + (orientation - 1) / 2);
                    obj.orientation_r(new_ind) = orientation;
                    obj.cell_ind_r(new_ind) = sub2ind(obj.grid.sz, face_index(1), face_index(2) + (orientation - 1) / 2);
                otherwise
                    error("Provided direction has to be a name of one of the directions, either 'x' or 'r'");
            end
        end
        
        function [coeffs_x, coeffs_r] = apply_boundary_condition_value(obj, coeffs_x, coeffs_r)
            % Heat conductivity at boundary (using assumption that its normal derivative at boundary is 0)
            k_x = obj.energy_component.k(obj.cell_ind_x);
            k_r = obj.energy_component.k(obj.cell_ind_r);

            % 1 / (h + k * 2 / dx)
            temp_coeff_x = 1 ./ (obj.heat_trans_coeff_x + k_x .* obj.dx_inv_times_2);
            temp_coeff_r = 1 ./ (obj.heat_trans_coeff_r + k_r .* obj.dr_inv_times_2);

            % X faces
            coeffs_x(obj.lin_ind_x_p1) = k_x .* obj.dx_inv_times_2 .* temp_coeff_x .* ((1 - obj.orientation_x) * 0.5);
            coeffs_x(obj.lin_ind_x_p2) = k_x .* obj.dx_inv_times_2 .* temp_coeff_x .* ((1 + obj.orientation_x) * 0.5);
            coeffs_x(obj.lin_ind_x_p3) = obj.heat_trans_coeff_x .* obj.environment_temp_x .* temp_coeff_x;

            % R faces
            coeffs_r(obj.lin_ind_r_p1) = k_r .* obj.dr_inv_times_2 .* temp_coeff_r .* ((1 - obj.orientation_r) * 0.5);
            coeffs_r(obj.lin_ind_r_p2) = k_r .* obj.dr_inv_times_2 .* temp_coeff_r .* ((1 + obj.orientation_r) * 0.5);
            coeffs_r(obj.lin_ind_r_p3) = obj.heat_trans_coeff_r .* obj.environment_temp_r .* temp_coeff_r;
        end

        function [coeffs_x, coeffs_r] = apply_boundary_condition_normal_derivative(obj, coeffs_x, coeffs_r)
            % Heat conductivity at boundary (using assumption that its normal derivative at boundary is 0)
            k_x = obj.energy_component.k(obj.cell_ind_x);
            k_r = obj.energy_component.k(obj.cell_ind_r);

            % orientation * h * 2 / dx / (h + k * 2 / dx)
            temp_coeff_x = obj.orientation_x .* obj.heat_trans_coeff_x .* obj.dx_inv_times_2 ./ (obj.heat_trans_coeff_x + k_x .* obj.dx_inv_times_2);
            temp_coeff_r = obj.orientation_r .* obj.heat_trans_coeff_r .* obj.dr_inv_times_2 ./ (obj.heat_trans_coeff_r + k_r .* obj.dr_inv_times_2);

            % X face
            coeffs_x(obj.lin_ind_x_p1) = temp_coeff_x .* ((1 - obj.orientation_x) * 0.5);
            coeffs_x(obj.lin_ind_x_p2) = temp_coeff_x .* ((1 + obj.orientation_x) * 0.5);
            coeffs_x(obj.lin_ind_x_p3) = -temp_coeff_x .* obj.environment_temp_x;

            % R faces
            coeffs_r(obj.lin_ind_r_p1) = temp_coeff_r .* ((1 - obj.orientation_r) * 0.5);
            coeffs_r(obj.lin_ind_r_p2) = temp_coeff_r .* ((1 + obj.orientation_r) * 0.5);
            coeffs_r(obj.lin_ind_r_p3) = -temp_coeff_r .* obj.environment_temp_r;
        end
    
        function converted_bd = convert_into_correction_boundary(obj)
            % NOT YET IMPLEMENTED
            converted_bd = obj.copy();
        end
    end
end