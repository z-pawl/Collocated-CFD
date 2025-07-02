classdef FixedValueBoundary < IBoundaryType
    properties
        grid Grid2D {mustBeScalarOrEmpty}
        value (1,1) double

        % Each row consists of: x-index, r-index, orientation
        boundary_faces_x (:,3) double
        boundary_faces_r (:,3) double
    end
    methods
        function obj = FixedValueBoundary(name, grid, value)
            arguments
                name (1,1) string
                grid (1,1) Grid2D
                value (1,1) double
            end
            obj.name = name;
            obj.grid = grid;
            obj.value = value;
        end

        function add_boundary_face(obj, direction, face_index, orientation)
            arguments
                obj
                direction (1,1) string
                face_index (1,2) double
                orientation (1,1) double
            end

            if orientation ~= -1 && orientation ~= 1
                error("The orientation has to be either -1 or 1");
            end

            switch direction
                case "x"
                    if (1 > face_index(1) || face_index(1) > obj.grid.sz(1)+1) || (1 > face_index(2) || face_index(2) > obj.grid.sz(2))
                        error("Provided face index is out of bounds");
                    end
                    obj.boundary_faces_x(size(obj.boundary_faces_x,1)+1,:) = [face_index orientation];
                case "r"
                    if (1 > face_index(1) || face_index(1) > obj.grid.sz(1)) || (1 > face_index(2) || face_index(2) > obj.grid.sz(2)+1)
                        error("Provided face index is out of bounds");
                    end
                    obj.boundary_faces_r(size(obj.boundary_faces_r,1)+1,:) = [face_index orientation];
                otherwise
                    error("Provided direction has to be a name of one of the directions, either 'x' or 'r'");
            end
        end
        
        function [coeffs_x, coeffs_r] = apply_boundary_condition_value(obj, coeffs_x, coeffs_r)
            % X faces
            faces_x = obj.boundary_faces_x;
            for k = 1:size(faces_x,1)
                coeffs_x(faces_x(k,1),faces_x(k,2),:) = [0 0 obj.value];
            end

            % R faces
            faces_r = obj.boundary_faces_r;
            for k = 1:size(faces_r,1)
                coeffs_r(faces_r(k,1),faces_r(k,2),:) = [0 0 obj.value];
            end
        end

        function [coeffs_x, coeffs_r] = apply_boundary_condition_normal_derivative(obj, coeffs_x, coeffs_r)
            % X faces
            faces_x = obj.boundary_faces_x;
            dx = obj.grid.dx(:,1);

            for k = 1:size(faces_x,1)
                i = faces_x(k,1);
                j = faces_x(k,2);
                orientation = faces_x(k,3);

                if orientation == 1
                    dx_inv = 2 / dx(i);
                    coeffs_x(i,j,:) = [0, dx_inv, -obj.value * dx_inv];
                else
                    dx_inv = -2 / dx(i-1);
                    coeffs_x(i,j,:) = [dx_inv, 0, -obj.value * dx_inv];
                end
            end

            % R faces
            faces_r = obj.boundary_faces_r;
            dr = obj.grid.dr(1,:);

            for k = 1:size(faces_r,1)
                i = faces_r(k,1);
                j = faces_r(k,2);
                orientation = faces_r(k,3);

                if orientation == 1
                    dr_inv = 2 / dr(j);
                    coeffs_r(i,j,:) = [0, dr_inv, -obj.value * dr_inv];
                else
                    dr_inv = -2 / dr(j-1);
                    coeffs_r(i,j,:) = [dr_inv, 0, -obj.value * dr_inv];
                end
            end
        end
    
        function converted_bd = convert_into_correction_boundary(obj)
            converted_bd = obj.copy();
            converted_bd.value = 0;
        end
    end
end