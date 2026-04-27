function field = solve(coeff, field, num_of_iters, sweep_dirs, relax_factor)
    sz = size(field);
    sweep_dirs_len = numel(sweep_dirs);

    one_minus_relax = 1 - relax_factor;

    % The system of equation is solved iteratively num_of_iters number of times
    for s = 1:num_of_iters
        % TDMA algorithm is used to calculate the new values for each of the columns seperately
        switch sweep_dirs(mod(s-1,sweep_dirs_len)+1)
            case 1
                % i = 1
                rhs = coeff(:,1,6) + coeff(:,1,5) .* field(:,2);
                field(:,1) = one_minus_relax .* field(:,1) + ... 
                    relax_factor * tdma(-coeff(:,1,2), coeff(:,1,1), -coeff(:,1,3), rhs);
                % Interior
                for i = 2:sz(2)-1
                    % Calculating the right hand side of the equation
                    rhs=coeff(:,i,6) + coeff(:,i,4) .* field(:,i-1) + coeff(:,i,5) .* field(:,i+1);
                    field(:,i) = one_minus_relax .* field(:,i) + ... 
                        relax_factor .* tdma(-coeff(:,i,2), coeff(:,i,1), -coeff(:,i,3), rhs);
                end
                % i = sz(2)
                rhs = coeff(:,sz(2),6) + coeff(:,sz(2),4) .* field(:,sz(2)-1);
                field(:,sz(2)) = one_minus_relax .* field(:,sz(2)) + ... 
                    relax_factor * tdma(-coeff(:,sz(2),2), coeff(:,sz(2),1), -coeff(:,sz(2),3), rhs);
            case 2
                % i = sz(2)
                rhs = coeff(:,sz(2),6) + coeff(:,sz(2),4) .* field(:,sz(2)-1);
                field(:,sz(2)) = one_minus_relax .* field(:,sz(2)) + ... 
                    relax_factor * tdma(-coeff(:,sz(2),2), coeff(:,sz(2),1), -coeff(:,sz(2),3), rhs);
                % Interior
                for i = sz(2)-1:(-1):2
                    % Calculating the right hand side of the equation
                    rhs=coeff(:,i,6) + coeff(:,i,4) .* field(:,i-1) + coeff(:,i,5) .* field(:,i+1);
                    field(:,i) = one_minus_relax .* field(:,i) + ... 
                        relax_factor .* tdma(-coeff(:,i,2), coeff(:,i,1), -coeff(:,i,3), rhs);
                end
                % i = 1
                rhs = coeff(:,1,6) + coeff(:,1,5) .* field(:,2);
                field(:,1) = one_minus_relax .* field(:,1) + ... 
                    relax_factor * tdma(-coeff(:,1,2), coeff(:,1,1), -coeff(:,1,3), rhs);
            case 3
                % i = 1
                rhs = coeff(1,:,6) + coeff(1,:,3) .* field(2,:);
                field(1,:) = one_minus_relax * field(1,:) + ...
                    relax_factor * tdma(-coeff(1,:,4), coeff(1,:,1), -coeff(1,:,5), rhs);
                % Interior
                for i=2:sz(1)-1
                    % Calculating the right hand side of the equation
                    rhs=coeff(i,:,6) + coeff(i,:,2) .* field(i-1,:) + coeff(i,:,3) .* field(i+1,:);
                    field(i,:) = one_minus_relax * field(i,:) + ...
                        relax_factor * tdma(-coeff(i,:,4), coeff(i,:,1), -coeff(i,:,5), rhs);
                end
                % i = sz(1)
                rhs = coeff(sz(1),:,6) + coeff(sz(1),:,2) .* field(sz(1)-1,:);
                field(sz(1),:) = one_minus_relax * field(sz(1),:) + ...
                    relax_factor * tdma(-coeff(sz(1),:,4), coeff(sz(1),:,1), -coeff(sz(1),:,5), rhs);
            case 4
                % i = sz(1)
                rhs = coeff(sz(1),:,6) + coeff(sz(1),:,2) .* field(sz(1)-1,:);
                field(sz(1),:) = one_minus_relax * field(sz(1),:) + ...
                    relax_factor * tdma(-coeff(sz(1),:,4), coeff(sz(1),:,1), -coeff(sz(1),:,5), rhs);
                % Interior
                for i=2:sz(1)-1
                    % Calculating the right hand side of the equation
                    rhs=coeff(i,:,6) + coeff(i,:,2) .* field(i-1,:) + coeff(i,:,3) .* field(i+1,:);
                    field(i,:) = one_minus_relax * field(i,:) + ...
                        relax_factor * tdma(-coeff(i,:,4), coeff(i,:,1), -coeff(i,:,5), rhs);
                end
                % i = 1
                rhs = coeff(1,:,6) + coeff(1,:,3) .* field(2,:);
                field(1,:) = one_minus_relax * field(1,:) + ...
                    relax_factor * tdma(-coeff(1,:,4), coeff(1,:,1), -coeff(1,:,5), rhs);
        end
    end
end