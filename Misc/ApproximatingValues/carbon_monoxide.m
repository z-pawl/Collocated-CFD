global A B C D E F G;

A = [25.56759	35.15070];
B = [6.096130	1.300095];
C = [4.054656	-0.205921];
D = [-2.671301 0.013550];
E = [0.131021	-3.282780];
F = [-118.0089	-127.8375];
G = [227.3665 231.7120];

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
    if T <= 1300
        h0 = calculate_H0(T/1000, 1);
    else
        h0 = calculate_H0(T/1000, 2);
    end
end

function g0 = c_G0(T)
    if T <= 1300
        g0 = calculate_G0(T/1000, 1);
    else
        g0 = calculate_G0(T/1000, 2);
    end
end

% Fitting the polynomials
co_C0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_C0(x, 1), 0.7, 1.3) + integral(@(x) polyval(ei, x) .* calculate_C0(x, 2), 1.3, 1.4);
    co_C0 = co_C0 + c * ei;
end

co_H0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_H0(x, 1), 0.7, 1.3) + integral(@(x) polyval(ei, x) .* calculate_H0(x, 2), 1.3, 1.4);
    co_H0 = co_H0 + c * ei;
end

co_S0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_S0(x, 1), 0.7, 1.3) + integral(@(x) polyval(ei, x) .* calculate_S0(x, 2), 1.3, 1.4);
    co_S0 = co_S0 + c * ei;
end


co_G0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_G0(x, 1), 0.7, 1.3) + integral(@(x) polyval(ei, x) .* calculate_G0(x, 2), 1.3, 1.4);
    co_G0 = co_G0 + c * ei;
end



figure(1);
plot(700:5:1400, arrayfun(@(x) polyval(co_H0, x/1000), 700:5:1400));
hold on;
plot(700:5:1400, arrayfun(@(x) c_H0(x), 700:5:1400));

figure(2);
plot(700:5:1400, arrayfun(@(x) polyval(co_G0, x/1000), 700:5:1400));
hold on;
plot(700:5:1400, arrayfun(@(x) c_G0(x), 700:5:1400));

figure(3);
plot(700:5:1400, arrayfun(@(x) 100 * (polyval(co_H0, x/1000) - c_H0(x)) / c_H0(x), 700:5:1400));

figure(4);
plot(700:5:1400, arrayfun(@(x) 100 * (polyval(co_G0, x/1000) - c_G0(x)) / c_G0(x), 700:5:1400));