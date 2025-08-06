function entropy_generation = entropy_generation_heat_transfer(grid, energy_component)
    arguments
        grid (1,1) Grid2D
        energy_component EnergyComponent
    end

    % Partial derivatives of temperature at faces [K/m]
    [coeffs_der_T_x, coeffs_der_T_r] = central_differencing_scheme(grid);
    [coeffs_der_T_x, coeffs_der_T_r] = energy_component.temp_bds.apply_boundary_condition_normal_derivative(coeffs_der_T_x, coeffs_der_T_r);
    [der_T_x, der_T_r] = evaluate_faces(coeffs_der_T_x, coeffs_der_T_r, energy_component.temp);

    % Temperature gradient at cell centers [K/m]
    gradient_T_x = (der_T_x(1:end-1,:) + der_T_x(2:end,:)) / 2;
    gradient_T_r = (der_T_r(:,1:end-1) + der_T_r(:,2:end)) / 2;

    % S_gen_local = -q * grad(T) / T ^ 2 = k * grad(T) ^ 2 / T ^ 2 = k * grad(T)^2 / T^2
    entropy_generation = energy_component.k .* (gradient_T_x .^ 2 + gradient_T_r .^ 2) ./ (energy_component.temp .^ 2);
end