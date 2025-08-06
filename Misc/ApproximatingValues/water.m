global A B C D E F G;

A = 30.09200;
B = 6.832514;
C = 6.793435;
D = -2.534480;
E = 0.082139;
F = -250.8810;
G = 223.3967;

function H0 = calculate_H0(t)
    global A B C D E F G;
    H0 = A * t + B * t .^ 2 / 2 + C * t .^ 3 / 3 + D * t .^ 4 / 4 - E ./ t + F;
end

function S0 = calculate_S0(t)
    global A B C D E F G;
    S0 = A * log(t) + B * t + C * t .^ 2 / 2 + D * t .^ 3 / 3 - E ./ (2 * t .^ 2) + G;
end

function G0 = calculate_G0(t)
    G0 = calculate_H0(t) - t .* calculate_S0(t);
end

function h0 = c_H0(T)
    h0 = calculate_H0(T/1000);
end

function g0 = c_G0(T)
    g0 = calculate_G0(T/1000);
end

% Fitting the polynomials
h2o_H0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_H0(x), 0.7, 1.4);
    h2o_H0 = h2o_H0 + c * ei;
end

h2o_G0 = 0;
for i = 1:size(orthonormal_basis,1)
    ei = orthonormal_basis(i,:);
    c = integral(@(x) polyval(ei, x) .* calculate_G0(x), 0.7, 1.4);
    h2o_G0 = h2o_G0 + c * ei;
end



figure(1);
plot(700:5:1400, arrayfun(@(x) polyval(h2o_H0, x/1000), 700:5:1400));
hold on;
plot(700:5:1400, arrayfun(@(x) c_H0(x), 700:5:1400));

figure(2);
plot(700:1:1400, arrayfun(@(x) polyval(h2o_G0, x/1000), 700:1:1400));
hold on;
plot(700:1:1400, arrayfun(@(x) c_G0(x), 700:1:1400));

figure(3);
plot(700:5:1400, arrayfun(@(x) 100 * (polyval(h2o_H0, x/1000) - c_H0(x)) / c_H0(x), 700:5:1400));

figure(4);
plot(700:5:1400, arrayfun(@(x) 100 * (polyval(h2o_G0, x/1000) - c_G0(x)) / c_G0(x), 700:5:1400));