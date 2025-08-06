%% SMR

dH0_SMR = (co_H0 + 3 * h2_H0) - (ch4_H0 + h2o_H0);
dG0_SMR = (co_G0 + 3 * h2_G0) - (ch4_G0 + h2o_G0);

figure(1);
plot(700:5:1400, arrayfun(@(x) polyval(dH0_SMR, x/1000), 700:5:1400));
figure(2);
plot(700:5:1400, arrayfun(@(x) polyval(dG0_SMR, x/1000), 700:5:1400));

%% WGS

dH0_WGS = (co2_H0 + h2_H0) - (co_H0 + h2o_H0);
dG0_WGS = (co2_G0 + h2_G0) - (co_G0 + h2o_G0);

figure(3);
plot(700:5:1400, arrayfun(@(x) polyval(dH0_WGS, x/1000), 700:5:1400));
figure(4);
plot(700:5:1400, arrayfun(@(x) polyval(dG0_WGS, x/1000), 700:5:1400));

%% Effectiveness

K_eq_SMR = arrayfun(@(x) exp(-1000 * polyval(dG0_SMR, x/1000) / (8.314 * x)), 700:5:1400);
figure(5);
plot(700:5:1400, K_eq_SMR);

SC = 2;
conversion_effectiveness = arrayfun(@(x) (x*(SC+1)^2-(SC+1)*sqrt(x^2*(SC-1)^2+12*x*SC))/((2*x-6)*(SC+1)^2), K_eq_SMR) / (1/(1+SC)) * 100;
figure(6);
plot(700:5:1400, conversion_effectiveness);