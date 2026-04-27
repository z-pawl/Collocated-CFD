function flux_limiter = van_leer_flux_limiter(derivative_face, derivative_upwind_face)
    arguments
        derivative_face (:,:) double
        derivative_upwind_face (:,:) double 
    end

    % Writing in place
    derivative_upwind_face = derivative_upwind_face ./ derivative_face; % r
    derivative_upwind_face(derivative_upwind_face < 0) = 0;

    flux_limiter = 2 * derivative_upwind_face ./ (1 + derivative_upwind_face);
    flux_limiter(derivative_face == 0) = 0;

    % if derivative_face ~= 0
    %     r = derivative_upwind_face / derivative_face;
    %     if r >= 0
    %         flux_limiter = 2 * r / (1 + r);
    %     else
    %         flux_limiter = 0;
    %     end
    % else
    %     flux_limiter = 0;
    % end
end