classdef (Abstract) IBoundaryType < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    properties
        name (1,1) string
    end
    methods (Abstract)
        add_boundary_face(obj, direction, face_index, orientation)
        [coeffs_x, coeffs_r] = apply_boundary_condition_value(obj, coeffs_x, coeffs_r)
        [coeffs_x, coeffs_r] = apply_boundary_condition_normal_derivative(obj, coeffs_x, coeffs_r)
        
        % Boundary conditions can be generally written in the form: a*y+b*y'=c
        % If a field is assumed to be a sum of current value y* and a correction y# the equation becomes
        % a*y*+b*y*'+a*y#+b*y#'=c
        % Since the current field is also constrained by the boundary condition the above equations simplifies to:
        % a*y#_b*y#'=0
        % The below function creates another boundary that is the boundary condition for the correction field
        converted_bd = convert_into_correction_boundary(obj)
    end
end