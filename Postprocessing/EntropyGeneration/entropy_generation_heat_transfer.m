function entropy_generation = entropy_generation_heat_transfer(grid, energy_component)
    arguments
        grid (1,1) Grid2D
        energy_component EnergyComponent
    end

    % Partial derivatives at faces
    [coeffs_part_der_x, coeff_part_der_r] = central_differencing_scheme(grid, energy_component.temp_bds);
    [part_der_at_faces_x, part_der_at_faces_r] = evaluate_faces(coeffs_part_der_x, coeff_part_der_r, energy_component.temp);

    % Temperature gradient at cell centers
    temp_grad_x = grid.face_area_x(2:end,:) .* part_der_at_faces_x(2:end,:) - grid.face_area_x(1:end-1,:) .* part_der_at_faces_x(1:end-1,:);
    temp_grad_r = grid.face_area_r(:,2:end) .* part_der_at_faces_r(:,2:end) - grid.face_area_r(:,1:end-1) .* part_der_at_faces_r(:,1:end-1);

    clear coeffs_part_der_x coeff_part_der_r part_der_at_faces_x part_der_at_faces_r;

    % S_gen = q * grad(1/T) = -k * grad(T) * grad(1/T) = k * grad(T)^2 / T^2
    entropy_generation = energy_component.k .* (temp_grad_x .^ 2 + temp_grad_r .^ 2) ./ (energy_component.temp .^ 2);
end