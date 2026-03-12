% Tridiagonal matrix algorithm
%     [ ⋱ ⋱       ]
%     [ ⋱ ⋱ c     ]
% A = [   ⋱ b ⋱   ]
%     [     a ⋱ ⋱ ]
%     [       ⋱ ⋱ ]
% b,d - nx1 vectors
% a,c - nx1 vectors - script will ignore the first and the last element of the vector respectively
function d = tdma(a,b,c,d)

    % Dimension of the matrix A
    n = length(b);

    % Forward sweep
    for i = 2:n
        w = a(i)/b(i-1);
        b(i) = b(i) - w*c(i-1);
        d(i) = d(i) - w*d(i-1);
    end


    % Backward substitution
    d(n) = d(n) / b(n);
    for i = (n-1):(-1):1
        d(i) = (d(i)-c(i)*d(i+1))/b(i);
    end
end

