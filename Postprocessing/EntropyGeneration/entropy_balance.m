function [rel_diff, s_gen, s_inlet, s_outlet, s_wall, s_inlet_conv, s_inlet_heat, s_inlet_diff] = entropy_balance(grid, energy, flow, species_manager, ent_gen_total)
% Universal gas constant [J/(mol*K)]
R = 8.314472;

% Reference pressure [Pa]
p_ref = 10^5;

ch4 = species_manager.species("CH4");
h2o = species_manager.species("H2O");
h2 = species_manager.species("H2");
co = species_manager.species("CO");
co2 = species_manager.species("CO2");

%% Entropy flux due to convection

% Polynomials of molar standard entropy [J/(mol*K)] in terms of temperature [K]
ch4_S0 = [-1.94991128455512e-12, 1.09474452431210e-08, -3.72128688651781e-05, 0.121177457943402, 154.585965536217];
h2o_S0 = [-7.06746263727913e-12, 3.83544978245327e-08, -8.70478509773900e-05, 0.128550749082486, 159.950878075594];
h2_S0 = [-7.20662147200438e-12, 3.93139020864608e-08, -8.81350482881648e-05, 0.117342974696574, 104.901546332004];
co_S0 = [-5.78838038824146e-12, 3.17167703576275e-08, -7.41426857456568e-05, 0.109450636244988, 173.298568106114];
co2_S0 = [-5.17096259801385e-12, 3.17813205704063e-08, -8.56319140061735e-05, 0.150895904430302, 177.428899696023];

% Polynomials of molar standard Gibbs free energy (standard chemical
% potential) [J/mol] in terms of temperatury [K]
ch4_G0 = [-6.65102121866067e-10, 8.05690857694305e-06, -0.0560868114381124, -156.885601038810, -78661.9683940239];
h2o_G0 = [-2.13305097598383e-09, 1.34786888731391e-05, -0.0482923241253257, -168.063329978702, -243553.843335180];
h2_G0 = [-2.22000062580425e-09, 1.35068801121399e-05, -0.0423259033842999, -113.208443620762, -1288.26799630224];
co_G0 = [-1.86499124717225e-09, 1.21357290674449e-05, -0.0418229966928410, -179.838087525415, -111984.026777027];
co2_G0 = [-2.73479118706640e-09, 1.81980560408625e-05, -0.0653171446899215, -182.325399802711, -397246.421493856];

% Inlet molar fractions [-]
SC = 2;
inlet_x_ch4 = 1 / (1 + SC);
inlet_x_h2o = SC / (1 + SC);
inlet_x_h2 = 0;
inlet_x_co = 0;
inlet_x_co2 = 0;

% Outlet molar fractions [-] (values are collected from cell centers due to the
% zero flux bd condition)
outlet_x_ch4 = ch4.Y(end, :) .* species_manager.M_mix(end, :) / ch4.M;
outlet_x_h2o = h2o.Y(end, :) .* species_manager.M_mix(end, :) / h2o.M;
outlet_x_h2 = h2.Y(end, :) .* species_manager.M_mix(end, :) / h2.M;
outlet_x_co = co.Y(end, :) .* species_manager.M_mix(end, :) / co.M;
outlet_x_co2 = co2.Y(end, :) .* species_manager.M_mix(end, :) / co2.M;

% Inlet and outlet pressure [Pa]
[lin_int_x, lin_int_r] = linear_interpolation_scheme(grid);
[lin_int_x, lin_int_r] = flow.p_bds.apply_boundary_condition_value(lin_int_x, lin_int_r);
[p_x, ~] = evaluate_faces(lin_int_x, lin_int_r, flow.p);
inlet_p = p_x(1, :);
outlet_p = p_x(end, :);

% Inlet and outlet temperature [K]
[lin_int_x, lin_int_r] = linear_interpolation_scheme(grid);
[lin_int_x, lin_int_r] = energy.temp_bds.apply_boundary_condition_value(lin_int_x, lin_int_r);
[T_x, T_r] = evaluate_faces(lin_int_x, lin_int_r, energy.temp);
inlet_T = T_x(1, :);
outlet_T = T_x(end, :);

% Inlet and outlet density [kg/m^3]
[lin_int_x, lin_int_r] = linear_interpolation_scheme(grid);
[lin_int_x, lin_int_r] = grid.domain_boundary.apply_boundary_condition_value(lin_int_x, lin_int_r);
[rho_x, ~] = evaluate_faces(lin_int_x, lin_int_r, flow.rho);
inlet_rho = rho_x(1, :);
outlet_rho = rho_x(end, :);

% Inlet and outler molar fluxes [mol/(m^2s)]
g_to_kg = 1 / 1000;
inlet_n_flux = inlet_rho .* flow.vx_faces(1, :) ./ (species_manager.M_mix(1, :) * g_to_kg);
outlet_n_flux = outlet_rho .* flow.vx_faces(end, :) ./ (species_manager.M_mix(end, :) * g_to_kg);

% Inlet and outlet specific molar entropy of species [J/(mol*K)]
% s_i = s^0_i - R ln(p / p_0) - R ln(x_i)
inlet_s_ch4 = polyval(ch4_S0, inlet_T) - R * log(inlet_p / p_ref) - R * log(inlet_x_ch4);
inlet_s_h2o = polyval(h2o_S0, inlet_T) - R * log(inlet_p / p_ref) - R * log(inlet_x_h2o);
inlet_s_h2 = polyval(h2_S0, inlet_T) - R * log(inlet_p / p_ref) - R * log(inlet_x_h2);
inlet_s_co = polyval(co_S0, inlet_T) - R * log(inlet_p / p_ref) - R * log(inlet_x_co);
inlet_s_co2 = polyval(co2_S0, inlet_T) - R * log(inlet_p / p_ref) - R * log(inlet_x_co2);

outlet_s_ch4 = polyval(ch4_S0, outlet_T) - R * log(outlet_p / p_ref) - R * log(outlet_x_ch4);
outlet_s_h2o = polyval(h2o_S0, outlet_T) - R * log(outlet_p / p_ref) - R * log(outlet_x_h2o);
outlet_s_h2 = polyval(h2_S0, outlet_T) - R * log(outlet_p / p_ref) - R * log(outlet_x_h2);
outlet_s_co = polyval(co_S0, outlet_T) - R * log(outlet_p / p_ref) - R * log(outlet_x_co);
outlet_s_co2 = polyval(co2_S0, outlet_T) - R * log(outlet_p / p_ref) - R * log(outlet_x_co2);

% Inlet and outlet specific molar entropy [J/(mol*K)]
% s = x_i * s_i
inlet_s = inlet_s_ch4 .* inlet_x_ch4 + inlet_s_h2o .* inlet_x_h2o;
outlet_s = outlet_s_ch4 .* outlet_x_ch4 + outlet_s_h2o .* outlet_x_h2o + outlet_s_h2 .* outlet_x_h2 + outlet_s_co .* outlet_x_co + outlet_s_co2 .* outlet_x_co2;

% Inlet and outlet entropy convection flux [J/(K*m^2*s)]
inlet_conv_s = inlet_s .* inlet_n_flux;
outlet_conv_s = outlet_s .* outlet_n_flux;

% Inlet and outlet entropy convection flux per face [W/K]
inlet_conv_s_face = inlet_conv_s .* grid.face_area_x(1, :);
outlet_conv_s_face = outlet_conv_s .* grid.face_area_x(end, :);

%% Entropy flux due to the heat conduction through wall and inlet
% Temperature gradient [K/m]
[cds_x, cds_r] = central_differencing_scheme(grid);
[cds_x, cds_r] = energy.temp_bds.apply_boundary_condition_normal_derivative(cds_x, cds_r);
[T_der_x, T_der_r] = evaluate_faces(cds_x, cds_r, energy.temp);

% Heat flux through wall and inlet [W/m^2]
q_r = -energy.k(:, end) .* T_der_r(:, end);
q_x = -energy.k(1, :) .* T_der_x(1, :);

% Entropy flux through wall and inlet [W/(K*m^2)]
wall_s = -q_r ./ T_r(:, end);
inlet_heat_s = q_x  ./ T_x(1, :);

% Entropy flux through wall and inlet per face [W/K]
wall_s_face = wall_s .* grid.face_area_r(:, end);
inlet_heat_s_face = inlet_heat_s .* grid.face_area_x(1, :);

%% Entropy flux due to the species diffusion through inlet
% Inlet chemical potential [J/mol]
inlet_chem_pot_ch4 = polyval(ch4_G0, inlet_T) + R * inlet_T .* log(inlet_p / p_ref);
inlet_chem_pot_h2o = polyval(h2o_G0, inlet_T) + R * inlet_T .* log(inlet_p / p_ref);
inlet_chem_pot_h2 = polyval(h2_G0, inlet_T) + R * inlet_T .* log(inlet_p / p_ref);
inlet_chem_pot_co = polyval(co_G0, inlet_T) + R * inlet_T .* log(inlet_p / p_ref);
inlet_chem_pot_co2 = polyval(co2_G0, inlet_T) + R * inlet_T .* log(inlet_p / p_ref);

% Inlet species molar diffusion fluxes [mol/s]
inlet_diff_flux_ch4 = (ch4.Jx(1,:) - ch4.Y(1,:) .* species_manager.Jx_sum(1,:)) / (ch4.M * g_to_kg);
inlet_diff_flux_h2o = (h2o.Jx(1,:) - h2o.Y(1,:) .* species_manager.Jx_sum(1,:)) / (h2o.M * g_to_kg);
inlet_diff_flux_h2 = (h2.Jx(1,:) - h2.Y(1,:) .* species_manager.Jx_sum(1,:)) / (h2.M * g_to_kg);
inlet_diff_flux_co = (co.Jx(1,:) - co.Y(1,:) .* species_manager.Jx_sum(1,:)) / (co.M * g_to_kg);
inlet_diff_flux_co2 = (co2.Jx(1,:) - co2.Y(1,:) .* species_manager.Jx_sum(1,:)) / (co2.M * g_to_kg);

% Entropy flux through inlet per face [W/K]
inlet_diff_s_face = - (inlet_chem_pot_ch4 .* inlet_diff_flux_ch4 + ...
    inlet_chem_pot_h2o .* inlet_diff_flux_h2o + ...
    inlet_chem_pot_h2 .* inlet_diff_flux_h2 + ...
    inlet_chem_pot_co .* inlet_diff_flux_co + ...
    inlet_chem_pot_co2 .* inlet_diff_flux_co2) ./ inlet_T;

%% Entropy generation and fluxes [W/K]
s_inlet_heat = sum(sum(inlet_heat_s_face));
s_inlet_diff = sum(sum(inlet_diff_s_face));
s_inlet_conv = sum(sum(inlet_conv_s_face));

s_inlet = s_inlet_conv + s_inlet_heat + s_inlet_diff;
s_outlet = sum(sum(outlet_conv_s_face));
s_wall = sum(sum(wall_s_face));

s_gen = s_outlet - s_inlet - s_wall;

%% Relative difference
rel_diff = (s_gen - ent_gen_total) / ent_gen_total;

