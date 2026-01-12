%% SMR

dH0_SMR = (co_H0 + 3 * h2_H0) - (ch4_H0 + h2o_H0);
dG0_SMR = (co_G0 + 3 * h2_G0) - (ch4_G0 + h2o_G0);

%% WGS

dH0_WGS = (co2_H0 + h2_H0) - (co_H0 + h2o_H0);
dG0_WGS = (co2_G0 + h2_G0) - (co_G0 + h2o_G0);
