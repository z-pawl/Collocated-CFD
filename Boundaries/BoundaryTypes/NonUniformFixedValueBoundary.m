classdef NonUniformFixedValueBoundary < IBoundaryType
    properties
        grid Grid2D {mustBeScalarOrEmpty}

        % Linearized indices
        lin_ind_x_p1 (:,1) double % x face, page 1 (i,j,1)
        lin_ind_x_p2 (:,1) double % x face, page 2 (i,j,2)
        lin_ind_x_p3 (:,1) double % x face, page 3 (i,j,3)

        lin_ind_r_p1 (:,1) double % r face, page 1 (i,j,1)
        lin_ind_r_p2 (:,1) double % r face, page 2 (i,j,2)
        lin_ind_r_p3 (:,1) double % r face, page 3 (i,j,3)

        % Coefficient values
        % Boundary value
        bd_val_x_p1 (:,1) double % x face, page 1 (i,j,1)
        bd_val_x_p2 (:,1) double % x face, page 2 (i,j,2)
        bd_val_x_p3 (:,1) double % x face, page 3 (i,j,3)

        bd_val_r_p1 (:,1) double % r face, page 1 (i,j,1)
        bd_val_r_p2 (:,1) double % r face, page 2 (i,j,2)
        bd_val_r_p3 (:,1) double % r face, page 3 (i,j,3)

        % Boundary normal derivative
        bd_der_x_p1 (:,1) double % x face, page 1 (i,j,1)
        bd_der_x_p2 (:,1) double % x face, page 2 (i,j,2)
        bd_der_x_p3 (:,1) double % x face, page 3 (i,j,3)

        bd_der_r_p1 (:,1) double % r face, page 1 (i,j,1)
        bd_der_r_p2 (:,1) double % r face, page 2 (i,j,2)
        bd_der_r_p3 (:,1) double % r face, page 3 (i,j,3)
    end
    methods
        function obj = NonUniformFixedValueBoundary(name, grid)
            arguments
                name (1,1) string
                grid (1,1) Grid2D
            end
            obj.name = name;
            obj.grid = grid;
        end

        function add_boundary_face(obj, direction, face_index, orientation, value)
            arguments
                obj
                direction (1,1) string
                face_index (1,2) double
                orientation (1,1) double
                value (1,1) double
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

                    % Boundary value coefficients
                    obj.bd_val_x_p1(new_ind) = 0;
                    obj.bd_val_x_p2(new_ind) = 0;
                    obj.bd_val_x_p3(new_ind) = value;

                    % Boundary normal derivative coefficients
                    dx = obj.grid.dx(:,1);
                    i = face_index(1);
                    if orientation == 1
                        dx_inv = 2 / dx(i);
                        obj.bd_der_x_p1(new_ind) = 0;
                        obj.bd_der_x_p2(new_ind) = dx_inv;
                        obj.bd_der_x_p3(new_ind) = -dx_inv * value;
                    else
                        dx_inv = -2 / dx(i-1);
                        obj.bd_der_x_p1(new_ind) = dx_inv;
                        obj.bd_der_x_p2(new_ind) = 0;
                        obj.bd_der_x_p3(new_ind) = -dx_inv * value;
                    end
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

                    % Boundary value coefficients
                    obj.bd_val_r_p1(new_ind) = 0;
                    obj.bd_val_r_p2(new_ind) = 0;
                    obj.bd_val_r_p3(new_ind) = value;

                    % Boundary normal derivative coefficients
                    dr = obj.grid.dr(1,:);
                    j = face_index(2);
                    if orientation == 1
                        dr_inv = 2 / dr(j);
                        obj.bd_der_r_p1(new_ind) = 0;
                        obj.bd_der_r_p2(new_ind) = dr_inv;
                        obj.bd_der_r_p3(new_ind) = -dr_inv * value;
                    else
                        dr_inv = -2 / dr(j-1);
                        obj.bd_der_r_p1(new_ind) = dr_inv;
                        obj.bd_der_r_p2(new_ind) = 0;
                        obj.bd_der_r_p3(new_ind) = -dr_inv * value;
                    end
                otherwise
                    error("Provided direction has to be a name of one of the directions, either 'x' or 'r'");
            end
        end
        
        function [coeffs_x, coeffs_r] = apply_boundary_condition_value(obj, coeffs_x, coeffs_r)
            % X faces
            coeffs_x(obj.lin_ind_x_p1) = obj.bd_val_x_p1;
            coeffs_x(obj.lin_ind_x_p2) = obj.bd_val_x_p2;
            coeffs_x(obj.lin_ind_x_p3) = obj.bd_val_x_p3;

            % R faces
            coeffs_r(obj.lin_ind_r_p1) = obj.bd_val_r_p1;
            coeffs_r(obj.lin_ind_r_p2) = obj.bd_val_r_p2;
            coeffs_r(obj.lin_ind_r_p3) = obj.bd_val_r_p3;
        end

        function [coeffs_x, coeffs_r] = apply_boundary_condition_normal_derivative(obj, coeffs_x, coeffs_r)
            % X faces
            coeffs_x(obj.lin_ind_x_p1) = obj.bd_der_x_p1;
            coeffs_x(obj.lin_ind_x_p2) = obj.bd_der_x_p2;
            coeffs_x(obj.lin_ind_x_p3) = obj.bd_der_x_p3;

            % R faces
            coeffs_r(obj.lin_ind_r_p1) = obj.bd_der_r_p1;
            coeffs_r(obj.lin_ind_r_p2) = obj.bd_der_r_p2;
            coeffs_r(obj.lin_ind_r_p3) = obj.bd_der_r_p3;
        end
    
        function converted_bd = convert_into_correction_boundary(obj)
            converted_bd = obj.copy();
            converted_bd.bd_val_x_p3(:) = 0;
            converted_bd.bd_val_r_p3(:) = 0;
            converted_bd.bd_der_x_p3(:) = 0;
            converted_bd.bd_der_r_p3(:) = 0;
        end
    end
end