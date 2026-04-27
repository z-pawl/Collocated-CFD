function [coeffs_x, coeffs_r] = TVD_scheme(grid, vel_x_faces, vel_r_faces, flux_limiter, field_face_derivatives_x, field_face_derivatives_r)
    arguments
        grid (1,1) Grid2D
        vel_x_faces (:,:) double
        vel_r_faces (:,:) double
        flux_limiter (1,1) function_handle
        field_face_derivatives_x (:,:) double 
        field_face_derivatives_r (:,:) double
    end
    sz = grid.sz;

    coeffs_x = zeros([(sz + [1 0]) 3]);
    coeffs_r = zeros([(sz + [0 1]) 3]);

    %% Coefficients in the x direction
    mask_pos_vel_x = vel_x_faces(2:sz(1),:) >= 0;
    upwind_face_derivatives_x = field_face_derivatives_x(3:sz(1)+1,:);
    temp_field_x = field_face_derivatives_x(1:sz(1)-1,:);
    upwind_face_derivatives_x(mask_pos_vel_x) = temp_field_x(mask_pos_vel_x);

    % Limiter
    temp_field_x = flux_limiter(field_face_derivatives_x(2:sz(1),:), upwind_face_derivatives_x);

    % Coefficient computations
    temp_field_x(mask_pos_vel_x) = 1 - temp_field_x(mask_pos_vel_x) / 2;
    temp_field_x(~mask_pos_vel_x) = temp_field_x(~mask_pos_vel_x) / 2;
    coeffs_x(2:sz(1),:,1) = temp_field_x;
    coeffs_x(2:sz(1),:,2) = 1 - temp_field_x;

    %% Coefficients in the r direction
    mask_pos_vel_r = vel_r_faces(:,2:sz(2)) >= 0;
    upwind_face_derivatives_r = field_face_derivatives_r(:,3:sz(2)+1);
    temp_field_r = field_face_derivatives_r(:,1:sz(2)-1);
    upwind_face_derivatives_r(mask_pos_vel_r) = temp_field_r(mask_pos_vel_r);

    % Limiter
    temp_field_r = flux_limiter(field_face_derivatives_r(:,2:sz(2)), upwind_face_derivatives_r);

    % Coefficient computations
    temp_field_r(mask_pos_vel_r) = 1 - temp_field_r(mask_pos_vel_r) / 2;
    temp_field_r(~mask_pos_vel_r) = temp_field_r(~mask_pos_vel_r) / 2;
    coeffs_r(:,2:sz(2),1) = temp_field_r;
    coeffs_r(:,2:sz(2),2) = 1 - temp_field_r;
end