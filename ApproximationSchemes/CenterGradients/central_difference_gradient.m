function [coeffs_x, coeffs_r] = central_difference_gradient(grid, field_bd_conditions)
    arguments
        grid (1,1) Grid2D
        field_bd_conditions (1,1) Boundaries
    end
    sz = grid.sz;

    coeffs_x = zeros([sz 6]);
    coeffs_r = zeros([sz 6]);

    % Inverse of the distance between cell centers along the x direction
    center_inv_distance_x = repmat(2 ./ ([0; grid.dx] + [grid.dx; 0]), 1, sz(2));
end