function [coeffs_x, coeffs_r] = upwind_scheme(grid, vel_x_faces, vel_r_faces)
    arguments
        grid (1,1) Grid2D
        vel_x_faces (:,:) double
        vel_r_faces (:,:) double
    end
    sz = grid.sz;

    coeffs_x = zeros([(sz + [1 0]) 3]);
    coeffs_r = zeros([(sz + [0 1]) 3]);

    coeffs_x(:,:,1) = double(vel_x_faces >= 0);
    coeffs_x(:,:,2) = double(vel_x_faces < 0);

    coeffs_r(:,:,1) = double(vel_r_faces >= 0);
    coeffs_r(:,:,2) = double(vel_r_faces < 0);
end