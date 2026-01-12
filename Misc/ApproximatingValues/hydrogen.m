global A B C D E F G;

A = [33.066178	18.563083];
B = [-11.363417	12.257357];
C = [11.432816	-2.859786];
D = [-2.772874 0.268238];
E = [-0.158558	1.977990];
F = [-9.980797	-1.147438];
G = [172.707974 156.288133];

function H0 = calculate_H0(t, i)
    global A B C D E F G;
    H0 = A(i) * t + B(i) * t .^2 / 2 + C(i) * t .^3 / 3 + D(i) * t .^ 4 / 4 - E(i) ./ t + F(i);
end

function S0 = calculate_S0(t, i)
    global A B C D E F G;
    S0 = A(i) * log(t) + B(i) * t + C(i) * t .^2 / 2 + D(i) * t .^ 3 / 3 - E(i) ./ (2 * t .^ 2) + G(i);
end

function G0 = calculate_G0(t,i)
    G0 = calculate_H0(t, i) - t .* calculate_S0(t, i);
end

function h0 = c_H0(T)
    if T <= 1000
        h0 = calculate_H0(T/1000, 1);
    else
        h0 = calculate_H0(T/1000, 2);
    end
end

function g0 = c_G0(T)
    if T <= 1000
        g0 = calculate_G0(T/1000, 1);
    else
        g0 = calculate_G0(T/1000, 2);
    end
end

% Fitting the polynomials
h2_H0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_H0(x, 1), 0.7, 1) + integral(@(x) polyval(ei, x) .* calculate_H0(x, 2), 1, 1.4);
    h2_H0 = h2_H0 + c * ei;
end

h2_G0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_G0(x, 1), 0.7, 1) + integral(@(x) polyval(ei, x) .* calculate_G0(x, 2), 1, 1.4);
    h2_G0 = h2_G0 + c * ei;
end