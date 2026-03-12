global A B C D E F G;

A = [24.99735	58.16639];
B = [55.18696	2.720074];
C = [-33.69137	-0.492289];
D = [7.948387 0.038844];
E = [-0.136638	-6.447293];
F = [-403.6075	-425.9186];
G = [228.2431 263.6125];

function C0 = calculate_C0(t, i)
    global A B C D E F G;
    C0 = A(i) + B(i) * t + C(i) * t .^ 2 + D(i) * t .^ 3 + E(i) ./ (t .^ 2);
end

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
    if T <= 1200
        h0 = calculate_H0(T/1000, 1);
    else
        h0 = calculate_H0(T/1000, 2);
    end
end

function g0 = c_G0(T)
    if T <= 1200
        g0 = calculate_G0(T/1000, 1);
    else
        g0 = calculate_G0(T/1000, 2);
    end
end

% Fitting the polynomials
co2_C0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_C0(x, 1), 0.7, 1.2) + integral(@(x) polyval(ei, x) .* calculate_C0(x, 2), 1.2, 1.4);
    co2_C0 = co2_C0 + c * ei;
end

co2_H0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_H0(x, 1), 0.7, 1.2) + integral(@(x) polyval(ei, x) .* calculate_H0(x, 2), 1.2, 1.4);
    co2_H0 = co2_H0 + c * ei;
end

co2_S0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_S0(x, 1), 0.7, 1.2) + integral(@(x) polyval(ei, x) .* calculate_S0(x, 2), 1.2, 1.4);
    co2_S0 = co2_S0 + c * ei;
end

co2_G0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_G0(x, 1), 0.7, 1.2) + integral(@(x) polyval(ei, x) .* calculate_G0(x, 2), 1.2, 1.4);
    co2_G0 = co2_G0 + c * ei;
end



figure(1);
plot(700:5:1400, arrayfun(@(x) polyval(co2_H0, x/1000), 700:5:1400));
hold on;
plot(700:5:1400, arrayfun(@(x) c_H0(x), 700:5:1400));

figure(2);
plot(700:5:1400, arrayfun(@(x) polyval(co2_G0, x/1000), 700:5:1400));
hold on;
plot(700:5:1400, arrayfun(@(x) c_G0(x), 700:5:1400));

figure(3);
plot(700:5:1400, arrayfun(@(x) 100 * (polyval(co2_H0, x/1000) - c_H0(x)) / c_H0(x), 700:5:1400));

figure(4);
plot(700:5:1400, arrayfun(@(x) 100 * (polyval(co2_G0, x/1000) - c_G0(x)) / c_G0(x), 700:5:1400));