function local_residual = calculate_local_residual(coeff, field)
    sz = size(field);

    % A column and a row of zeros used to pad the data in order to make sure that there is no index out of bounds error
    col0 = zeros(sz(1),1);
    row0 = zeros(1,sz(2));

    % r = aW*TW+aE*TE+aS*TS+aN*TN+b-aP*TP
    local_residual = [row0; coeff(2:end,:,2) .* field(1:end-1,:)] + [coeff(1:end-1,:,3) .* field(2:end,:); row0] ... 
        + [col0 coeff(:,2:end,4) .* field(:,1:end-1)] + [coeff(:,1:end-1,5) .* field(:,2:end) col0] ...
        + coeff(:,:,6) - coeff(:,:,1) .* field;
end