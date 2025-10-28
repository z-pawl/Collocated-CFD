function inner_prod = inner_product(v1, v2, p, q)
    v1_v2_int = polyint(conv(v1, v2));
    inner_prod = polyval(v1_v2_int, q) - polyval(v1_v2_int, p);
end

function vnorm = poly_norm(v, p, q)
    vnorm = sqrt(inner_product(v, v, p, q));
end

% Generates an orthonormal basis over a specified range (p:q) using the
% modified gram schmidt procedure
function basis = modified_gram_schmidt(basis, p, q)
    function orthonormal_basis = modified_gram_schmidt_helper(original_basis, p, q)
        n = size(original_basis, 1);
        orthonormal_basis = original_basis;
    
        for t = 1:n
            v = original_basis(t,:);
    
            for s = 1:t-1
                ej = orthonormal_basis(s,:);
    
                % Projection factor
                r = inner_product(v, ej, p, q) / inner_product(ej, ej, p, q);
    
                v = v - r * ej;
            end
    
            orthonormal_basis(t,:) = v / poly_norm(v, p, q);
        end
    end

    n = size(basis, 1);
    err = zeros(n) - eye(n);
    noi = 0;
    
    while sum(sum(abs(err))) >= 1e-7
        basis = modified_gram_schmidt_helper(basis, p, q);
        for i = 1:n
            for j = 1:n
                err(i, j) = inner_product(basis(i,:), basis(j,:), p, q);
            end
        end
        err = err - eye(n)
        noi = noi + 1;
    end

    fprintf("The calculation took %d iterations\n", noi);
end

n = 4;
p = 0.7;
q = 1.400;

original_basis = zeros(n+1);
for i = 1:n+1
    original_basis(i, n+2-i) = 1;
end

orthonormal_basis = modified_gram_schmidt(original_basis, p, q);

figure(1);
hold on;

for i = 0:n
    ei = orthonormal_basis(i+1, :);
    plot(p:0.005:q, arrayfun(@(x) polyval(ei, x), p:0.005:q), DisplayName = "e" + i);
end

clear ei i n original_basis p q;