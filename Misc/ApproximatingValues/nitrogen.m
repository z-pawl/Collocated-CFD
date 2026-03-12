global A B C D E F G;

A = 19.50583;
B = 19.88705;
C = -8.598535;
D = 1.369784;
E = 0.527601;
F = -4.935202;
G = 212.3900;

function C0 = calculate_C0(t)
    global A B C D E F G;
    C0 = A + B * t + C * t .^ 2 + D * t .^ 3 + E ./ (t .^ 2);
end

function H0 = calculate_H0(t)
    global A B C D E F G;
    H0 = A * t + B * t .^2 / 2 + C * t .^3 / 3 + D * t .^ 4 / 4 - E ./ t + F;
end

function S0 = calculate_S0(t)
    global A B C D E F G;
    S0 = A * log(t) + B * t + C * t .^2 / 2 + D * t .^ 3 / 3 - E ./ (2 * t .^ 2) + G;
end

function G0 = calculate_G0(t)
    G0 = calculate_H0(t) - t .* calculate_S0(t);
end

% Fitting the polynomials
n2_C0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_C0(x), 0.7, 1.4);
    n2_C0 = n2_C0 + c * ei;
end

n2_H0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_H0(x), 0.7, 1.4);
    n2_H0 = n2_H0 + c * ei;
end

n2_S0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_S0(x), 0.7, 1.4);
    n2_S0 = n2_S0 + c * ei;
end

n2_G0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_G0(x), 0.7, 1.4);
    n2_G0 = n2_G0 + c * ei;
end