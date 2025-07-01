function coeff = implicit_relaxation(coeff, field, relaxation_factor)
    coeff(:,:,6) = coeff(:,:,6) + (1 - relaxation_factor) / relaxation_factor * coeff(:,:,1) .* field;
    coeff(:,:,1) = coeff(:,:,1) / relaxation_factor;
end