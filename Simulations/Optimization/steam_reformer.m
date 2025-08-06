function [entropy_generation ch4_in, ch4_out, h2_out] = steam_reformer(temp_inlet, temp_wall, v_inlet, SC)

%% Parametry symulacji

% Geometria
% Długość reformera [m]
length = 0.4;
% Promień reformera [m]
radius = 0.05;


% Parametry reformera

% Porowatość [-]
porosity = 0.65;
% Przepuszczalność [m^2]
permeability = 1e-7;
% Współczynnik bezwładności [-]
inertia_coefficient = 0.088;
% Przewodność cieplna fazy stałej [W/(mK)]
solid_phase_thermal_conductivity = 20.0;
% Gęstość katalizatora [g/m^3]
catalyst_density_val = 2.5e5;
% Symuluje się również 5 cm przed i za reformerem w celu lepszego ukazania
% dyfuzji
catalyst_density = zeros(2 * [40 15]);
catalyst_density(11:70,:) = catalyst_density_val;



% Stałe reakcji chemicznych

% Rząd reakcji SMR względem metanu [-]
a = 0.89;
% Rząd reakcji SMR względem pary wodnej [-]
b = 0.05;
% Czynnik przedwykładniczy reakcji SMR [mol/(s*g*Pa^(a+b))]
A_st = 1.354e-3;
% Energia aktywacji reakcji SMR [J/mol]
E_a = 122500;
% Entalpia reakcji SMR [J/mol] (funkcja (wielomian) temperatury [K])
delta_H_st_f = [-1.515651567024407e-09, 1.538155721897593e-05, -0.055034658100429, 80.715671242009880, 1.861716447722712e+05];
% Entalpia reakcji WGS [J/mol] (funkcja (wielomian) temperatury [K])
delta_H_sh_f = [2.563977875923532e-09, -1.163445257819235e-05, 0.017949454968009, -1.182340818199378, -4.247549893130190e+04];


% Warunki brzegowe

% Temperatura na wlocie [K]
temp_inlet = 950;
% Temperatura ścian [K]
temp_wall = 1100;
% Prędkość na wlocie [m/s]
v_inlet = 0.25;
% Ciśnienie atmosferyczne (ciśnienie na wylocie) [Pa]
p_atm = 101325;
% Steam-carbon ratio na wlocie [-]
SC = 2;


% Parametry gazów

% Stała, której użycie wynika z użycia we wzorach na parametry jako
% argumentu temperatury podzielonej przez tysiąc
K = 1e-3;
% Stała w celu konwersji jednostek z uP na Pa*s
uP_to_Pas = 1e-7;

% Każdy z gazów jest opisany poprzez następujące wartości:
% Masa molowa [g/mol]
% Przewodność cieplna czystej substancji - opisana jako wielomian temperatury w [K]; [J/(mK)]
% Lepkość dynamiczna czystej substancji - opisana jako wielomian temperatury w [K]; [Pa*s]
% Molowa pojemność cieplna czystej substancji - opisana jako wielomian temperatury w [K]; [J/(mol*K)]
% Objętość dyfuzyjna według Fuller et al. [?]

% Wodór (H2)
M_H2 = 2.016;
therm_cond_H2 = 0.01 * [-1.8954*K^6 11.972*K^5 -31.939*K^4 47.763*K^3 -47.190*K^2 62.892*K 1.5040];
visc_H2 = uP_to_Pas*[-9.9892*K^6 62.966*K^5 -167.51*K^4 249.41*K^3 -244.34*K^2 299.78*K 15.553];
heat_cap_H2 = [-6.4725*K^6 46.903*K^5 -136.15*K^4 199.29*K^3 -150.55*K^2 56.036*K 21.157];
diffusion_volume_H2 = 6.12;

% Tlenek węgla (CO)
M_CO = 28.010;
therm_cond_CO = 0.01 * [-2.3224*K^6 13.379*K^5 -30.818*K^4 36.018*K^3 -23.186*K^2 13.999*K -0.2815];
visc_CO = uP_to_Pas*[-32.298*K^6 208.42*K^5 -572.14*K^4 883.75*K^3 875.90*K^2 793.65*K -4.9137];
heat_cap_CO = [-7.6538*K^6 37.756*K^5 -66.346*K^4 41.974*K^3 5.2062*K^2 -8.1781*K 30.429];
diffusion_volume_CO = 18.0;

% Dwutlenek węgla (CO2)
M_CO2 = 44.009;
therm_cond_CO2 = 0.01 * [18.698*K^6 -101.12*K^5 216.83*K^4 -233.29*K^3 129.65*K^2 -27.018*K 2.8888];
visc_CO2 = uP_to_Pas*[-0.4564*K^6 14.450*K^5 -85.929*K^4 244.22*K^3 -432.49*K^2 680.07*K -20.434];
heat_cap_CO2 = [-35.992*K^6 214.58*K^5 -519.9*K^4 657.88*K^3 -471.33*K^2 204.60*K 4.3669];
diffusion_volume_CO2 = 26.7;

% Metan (CH4)
M_CH4 = 16.043;
therm_cond_CH4 = 0.01 * [3.2774*K^6 -17.283*K^5 38.251*K^4 -47.440*K^3 37.413*K^2 1.8732*K 0.4796];
visc_CH4 = uP_to_Pas*[-22.920*K^6 140.48*K^5 -367.06*K^4 548.11*K^3 -543.82*K^2 529.37*K -9.9989];
heat_cap_CH4 = [61.321*K^6 -358.75*K^5 856.93*K^4 -1068.7*K^3 712.55*K^2 -178.59*K 47.964];
diffusion_volume_CH4 = 25.14;

% Para wodna (H2O)
M_H2O = 18.015;
therm_cond_H2O = 0.01 * [4.1531*K^6 -18.974*K^5 35.993*K^4 -41.390*K^3 35.922*K^2 -7.9139*K 2.0103];
visc_H2O = uP_to_Pas*[19.591*K^6 -126.96*K^5 348.12*K^4 -522.38*K^3 419.50*K^2 244.93*K -6.7541];
heat_cap_H2O = [14.015*K^6 -79.409*K^5 181.54*K^4 -217.08*K^3 146.01*K^2 -41.205*K 37.373];
diffusion_volume_H2O = 13.1;



%% Parametry solvera

% Ogólne

% Maksymalna liczba iteracji
max_iters = 10000;
% Rozmiar siatki
sz = 2 * [40 15];
% Rozmiary elementów
dx = length * ones(sz(1), 1) / sz(1);
dr = radius * ones(1, sz(2)) / sz(2);


% Przepływ

% Liczba wewnętrznych iteracji na iterację zewnętrzną dla równań przepływu
inner_iters_flow = 3;
% Liczba iteracji solvera na wewnętrzną iterację dla równań Navier'a-Stokes'a
solver_iters_v = 1;
% Liczba iteracji solvera na wewnętrzną iterację dla równania poprawki ciśnienia
solver_iters_p = 70;
% Kryterium konwergencji dla równań Navier'a-Stokes'a
tol_v = 3e-6;
% Kryterium konwergencji dla ciśnienia
tol_p = 5e-6;
% Współczynnik relaksacji prędkości
rel_fact_v = 0.75;
% Współczynnik relaksacji ciśnienia
rel_fact_p = 0.25;


% Energia

% Liczba wewnętrznych iteracji na iterację zewnętrzną dla równania energii
inner_iters_energy = 1;
% Liczba iteracji solvera na wewnętrzną iterację dla równania energii
solver_iters_energy = 15;
% Kryterium konwergencji dla równania energii
tol_energy = 1e-7;
% Współczynnik relaksacji temperatury
rel_fact_temp = 0.8;


% Substancje chemiczne

% Liczba wewnętrznych iteracji na iterację zewnętrzną dla substancji chemicznych
inner_iters_species = 3;
% Liczba iteracji solvera na wewnętrzną iterację dla każdej substancji chemicznej
solver_iters_species = 30;
% Kryterium konwergencji dla każdej substancji chemicznej
tol_species = 1e-7;
% Współczynnik relaksacji tempa reakcji SMR
rel_fact_R_st = 0.006;
% Współczynnik relaksacji tempa reakcji WGS
rel_fact_R_sh = 0.006;
% Współczynnik relaksacji ułamków masowych dla każdej substancji chemicznej
rel_fact_Y = 0.8;



%% Definicja komponentów

% Inicjalizacja komponentów

% Siatka
offset = 1e-30; % Niewielki offset, aby zapobiec dzieleniu przez zero
grid = Grid2D(dx, dr, offset);

% Komponent przepływu
flow = FlowComponent(grid, v_inlet * ones(sz), zeros(sz), p_atm * ones(sz), porosity * ones(sz), @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), tol_v, tol_p, rel_fact_v, rel_fact_p, inner_iters_flow, solver_iters_v, solver_iters_p);

% Komponent energii
energy = EnergyComponent(grid, flow, temp_inlet * ones(sz), @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), tol_energy, rel_fact_temp, inner_iters_energy, solver_iters_energy);

% Substancje chemiczne
% Menadżer
species_manager = SpeciesManagerComponent(grid, flow, energy, @(x) 0, a, b, A_st, E_a, catalyst_density, rel_fact_R_st, rel_fact_R_sh, inner_iters_species);
% Poszczególne substancje
% Wodór (H2)
h2 = SpeciesComponent("H2", grid, flow, species_manager, zeros(sz), heat_cap_H2, visc_H2, therm_cond_H2, M_H2, diffusion_volume_H2, tol_species, rel_fact_Y, solver_iters_species);
% Tlenek węgla (CO)
co = SpeciesComponent("CO", grid, flow, species_manager, zeros(sz), heat_cap_CO, visc_CO, therm_cond_CO, M_CO, diffusion_volume_CO, tol_species, rel_fact_Y, solver_iters_species);
% Dwutlenek węgla (CO2)
co2 = SpeciesComponent("CO2", grid, flow, species_manager, zeros(sz), heat_cap_CO2, visc_CO2, therm_cond_CO2, M_CO2, diffusion_volume_CO2, tol_species, rel_fact_Y, solver_iters_species);
% Metan (CH4)
ch4 = SpeciesComponent("CH4", grid, flow, species_manager, M_CH4/(M_CH4+SC*M_H2O) * ones(sz), heat_cap_CH4, visc_CH4, therm_cond_CH4, M_CH4, diffusion_volume_CH4, tol_species, rel_fact_Y, solver_iters_species);
% Para wodna (H2O)
h2o = SpeciesComponent("H2O", grid, flow, species_manager, SC*M_H2O/(M_CH4+SC*M_H2O) * ones(sz), heat_cap_H2O, visc_H2O, therm_cond_H2O, M_H2O, diffusion_volume_H2O, tol_species, rel_fact_Y, solver_iters_species);

% Menadżer właściwości termofizycznych
properties_manager = ThermophysicalProperties(grid, flow, energy, species_manager, solid_phase_thermal_conductivity);


% Właściwości termofizyczne
% Komponent przepływu
flow.rho_function = @properties_manager.rho_function;
flow.visc_function = @properties_manager.visc_function;
% Komponent energii
energy.k_function = @properties_manager.k_function;
energy.cp_function = @properties_manager.cp_function;
% Komponent substancji chemicznych
species_manager.D_f = @properties_manager.D_function;


% Człony źródłowe
% Źródła komponentu prędkości
% Liniowy człon źródłowy dla prędkości w kierunku osiowym
function src_lin_vx = lin_src_vx(flow, permeability, inertia_coefficient)
    % Człon źródłowy = vx*(-μ/Kp-rho*f/sqrt(Kp)*sqrt(vx^2+vr^2))
    src_lin_vx = -flow.visc / permeability - flow.rho * inertia_coefficient / sqrt(permeability) .* sqrt(flow.vx .^ 2 + flow.vr .^ 2);
end
% Liniowy człon źródlowy dla prędkości w kierunku promieniowym
function src_lin_vr = lin_src_vr(flow, permeability, inertia_coefficient)
    % Człon źródłowy = vr*(-μ/Kp-rho*f/sqrt(Kp)*sqrt(vx^2+vr^2))
    src_lin_vr = -flow.visc / permeability - flow.rho * inertia_coefficient / sqrt(permeability) .* sqrt(flow.vx .^ 2 + flow.vr .^ 2);
end

flow.src_linx_function = @() lin_src_vx(flow, permeability, inertia_coefficient);
flow.src_linr_function = @() lin_src_vr(flow, permeability, inertia_coefficient);

% Źródła komponentu energii
function src_q = q_source(energy, species_manager, delta_H_st_f, delta_H_sh_f)
    delta_H_st = polyval(delta_H_st_f, energy.temp);
    delta_H_sh = polyval(delta_H_sh_f, energy.temp);
    src_q = -delta_H_st .* species_manager.R_st - delta_H_sh .* species_manager.R_sh;
end

energy.q_function = @() q_source(energy, species_manager, delta_H_st_f, delta_H_sh_f);


% Warunki brzegowe
% Wlot
inlet_vx = FixedValueBoundary("Inlet vx", grid, v_inlet);
inlet_vr = FixedValueBoundary("Inlet vr", grid, 0);
% Warunek brzegowy ciśnienia został wyznaczony na podstawie prawa
% Darcy'ego-Forchheimer'a - grad(p) = -μ/K*v - f*ρ/sqrt(K)*v*|v|
% Należy wyznaczyć lepkość i gęstość na wlocie (efekt zmiany ciśnienia
% został pominięty)
species_manager.update_properties();
flow.update_properties();
visc_inlet = flow.visc(1,1);
rho_inlet = flow.rho(1,1);
p_der_inlet = -visc_inlet / permeability * v_inlet - inertia_coefficient * rho_inlet / sqrt(permeability) * v_inlet ^ 2;
inlet_p = FixedNormalDerivativeBoundary("Inlet p", grid, p_der_inlet);
clear visc_inlet rho_inlet p_der_inlet;
inlet_T = FixedValueBoundary("Inlet T", grid, temp_inlet);
inlet_ch4 = FixedValueBoundary("Inlet CH4", grid, M_CH4/(M_CH4+SC*M_H2O));
inlet_h2o = FixedValueBoundary("Inlet H2O", grid, SC*M_H2O/(M_CH4+SC*M_H2O));
inlet_rest_species = FixedValueBoundary("Inlet rest species", grid, 0);

% Wylot
outlet_vx = FixedNormalDerivativeBoundary("Outlet vx", grid, 0);
outlet_vr = FixedNormalDerivativeBoundary("Outlet vr", grid, 0);
outlet_p = FixedValueBoundary("Outlet p", grid, p_atm);
outlet_T = FixedNormalDerivativeBoundary("Outlet T", grid, 0);
outlet_species = FixedNormalDerivativeBoundary("Outlet species", grid, 0);

% Oś
axis_vx = FixedNormalDerivativeBoundary("Axis vx", grid, 0);
axis_vr = FixedValueBoundary("Axis vr", grid, 0);
axis_p = FixedNormalDerivativeBoundary("Axis p", grid, 0);
axis_T = FixedNormalDerivativeBoundary("Axis T", grid, 0);
axis_species = FixedNormalDerivativeBoundary("Axis species", grid, 0);

% Ściana
wall_vx = FixedValueBoundary("Wall vx", grid, 0);
wall_vr = FixedValueBoundary("Wall vr", grid, 0);
wall_p = FixedNormalDerivativeBoundary("Wall p", grid, 0);
wall_T = FixedValueBoundary("Wall T", grid, temp_wall);
wall_species = FixedNormalDerivativeBoundary("Wall species", grid, 0);

% Granica
domain_boundary = FixedNormalDerivativeBoundary("Domain boundary", grid, 0);

for j = 1:sz(2)
    % Wlot
    inlet_vx.add_boundary_face("x", [1 j], 1);
    inlet_vr.add_boundary_face("x", [1 j], 1);
    inlet_p.add_boundary_face("x", [1 j], 1);
    inlet_T.add_boundary_face("x", [1 j], 1);
    inlet_ch4.add_boundary_face("x", [1 j], 1);
    inlet_h2o.add_boundary_face("x", [1 j], 1);
    inlet_rest_species.add_boundary_face("x", [1 j], 1);
    domain_boundary.add_boundary_face("x", [1 j], 1);

    % Wylot
    outlet_vx.add_boundary_face("x", [sz(1)+1 j], -1);
    outlet_vr.add_boundary_face("x", [sz(1)+1 j], -1);
    outlet_p.add_boundary_face("x", [sz(1)+1 j], -1);
    outlet_T.add_boundary_face("x", [sz(1)+1 j], -1);
    outlet_species.add_boundary_face("x", [sz(1)+1 j], -1);
    domain_boundary.add_boundary_face("x", [sz(1)+1 j], -1);
end

for i = 1:sz(1)
    % Axis
    axis_vx.add_boundary_face("r", [i 1], 1);
    axis_vr.add_boundary_face("r", [i 1], 1);
    axis_p.add_boundary_face("r", [i 1], 1);
    axis_T.add_boundary_face("r", [i 1], 1);
    axis_species.add_boundary_face("r", [i 1], 1);
    domain_boundary.add_boundary_face("r", [i 1], 1);

    % Wall
    wall_vx.add_boundary_face("r", [i sz(2)+1], -1);
    wall_vr.add_boundary_face("r", [i sz(2)+1], -1);
    wall_p.add_boundary_face("r", [i sz(2)+1], -1);
    wall_T.add_boundary_face("r", [i sz(2)+1], -1);
    wall_species.add_boundary_face("r", [i sz(2)+1], -1);
    domain_boundary.add_boundary_face("r", [i sz(2)+1], -1);
end

grid.domain_boundary.add_boundary(domain_boundary);

flow.vx_bds.add_boundary(inlet_vx);
flow.vx_bds.add_boundary(outlet_vx);
flow.vx_bds.add_boundary(axis_vx);
flow.vx_bds.add_boundary(wall_vx);

flow.vr_bds.add_boundary(inlet_vr);
flow.vr_bds.add_boundary(outlet_vr);
flow.vr_bds.add_boundary(axis_vr);
flow.vr_bds.add_boundary(wall_vr);

flow.p_bds.add_boundary(inlet_p);
flow.p_bds.add_boundary(outlet_p);
flow.p_bds.add_boundary(axis_p);
flow.p_bds.add_boundary(wall_p);

flow.convert_pressure_bd_conditions();

energy.temp_bds.add_boundary(inlet_T);
energy.temp_bds.add_boundary(outlet_T);
energy.temp_bds.add_boundary(axis_T);
energy.temp_bds.add_boundary(wall_T);

ch4.species_bds.add_boundary(inlet_ch4);
ch4.species_bds.add_boundary(outlet_species);
ch4.species_bds.add_boundary(axis_species);
ch4.species_bds.add_boundary(wall_species);

h2o.species_bds.add_boundary(inlet_h2o);
h2o.species_bds.add_boundary(outlet_species);
h2o.species_bds.add_boundary(axis_species);
h2o.species_bds.add_boundary(wall_species);

co.species_bds.add_boundary(inlet_rest_species);
co.species_bds.add_boundary(outlet_species);
co.species_bds.add_boundary(axis_species);
co.species_bds.add_boundary(wall_species);

co2.species_bds.add_boundary(inlet_rest_species);
co2.species_bds.add_boundary(outlet_species);
co2.species_bds.add_boundary(axis_species);
co2.species_bds.add_boundary(wall_species);

h2.species_bds.add_boundary(inlet_rest_species);
h2.species_bds.add_boundary(outlet_species);
h2.species_bds.add_boundary(axis_species);
h2.species_bds.add_boundary(wall_species);


% Inicjalizacja
species_manager.update_properties();
energy.update_properties();
flow.update_properties();
[coeff_vx, coeff_vr] = flow.get_coefficients_v();
flow.update_face_velocities(flow.vx, flow.vr, coeff_vx(:,:,1), coeff_vr(:,:,1));
clear coeff_vx voeff_vr;



%% Symulacja
% Numer iteracji
noi = 0;
converged = false;
crashed = false;
species_names = keys(species_manager.species);
num_species = numel(species_names);
while ~converged && noi < max_iters
    noi = noi + 1;

    % Iterowanie każdego komponentu
    converged = true;
    converged = flow.iterate() && converged;
    converged = energy.iterate() && converged;
    converged = species_manager.iterate() && converged;

    if mod(noi, 500) == 0
        species_manager.urf_J_sum = clip(25 * species_manager.urf_J_sum, 0, 1);
    end

    % % Pokazywanie rezyduów
    % if mod(noi, 5) == 0
    %     fprintf('Current iteration: %d\n', noi);
    %     fprintf('The vx residual is equal to: %d\n', flow.residual_history_vx(noi));
    %     fprintf('The vr residual is equal to: %d\n', flow.residual_history_vr(noi));
    %     fprintf('The continuity residual is equal to: %d\n', flow.residual_history_continuity(noi));
    %     fprintf('The T residual is equal to: %d\n', energy.residual_history(noi));
    %     for i = 1:num_species
    %         res_hist = species_manager.residual_history(species_names(i));
    %         fprintf('The %s residual is equal to: %d\n', species_names(i), res_hist(noi));
    %     end
    % end

    % Jeśli jakieś pole posiada rezydua równe NaN lub Inf przez cały czas w
    % trakcie ostatnich 20 iteracji to zakłada się, że symulacja
    % zdywergowała
    if noi > 20
        crashed = crashed || all(isnan(flow.residual_history_vx(noi-20:noi)) | isinf(flow.residual_history_vx(noi-20:noi)));
        crashed = crashed || all(isnan(flow.residual_history_vr(noi-20:noi)) | isinf(flow.residual_history_vr(noi-20:noi)));
        crashed = crashed || all(isnan(flow.residual_history_continuity(noi-20:noi)) | isinf(flow.residual_history_continuity(noi-20:noi)));
        crashed = crashed || all(isnan(energy.residual_history(noi-20:noi)) | isinf(energy.residual_history(noi-20:noi)));
        for i = 1:num_species
            res_hist = species_manager.residual_history(species_names(i));
            crashed = crashed || all(isnan(res_hist(noi-20:noi)) | isinf(res_hist(noi-20:noi)));
        end
        
    end

    if crashed
        break;
    end
end



%% Przedstawienie wyników

% % Siatka
% [X, R] = meshgrid(grid.cent_pos_x, grid.cent_pos_r);
% 
% % Wykres wartości prędkości
% figure(1);
% colormap(jet);
% contourf(X, R, transpose(sqrt(flow.vx .^ 2 + flow.vr .^ 2)), 30);
% pbaspect([length radius 1]);
% 
% % Wykres prędkości w kierunku X
% figure(2);
% colormap(jet);
% contourf(X, R, transpose(flow.vx), 30);
% pbaspect([length radius 1]);
% 
% % Wykres prędkości w kierunku R
% figure(3);
% colormap(jet);
% contourf(X, R, transpose(flow.vr), 30);
% pbaspect([length radius 1]);
% 
% % Wykres ciśnienia
% figure(4);
% colormap(jet);
% contourf(X, R, transpose(flow.p), 30);
% pbaspect([length radius 1]);
% 
% % Wykres linii prądu
% figure(5);
% streamslice(X, R, transpose(flow.vx), transpose(flow.vr), 2);
% pbaspect([length radius 1]);
% 
% % Wykres temperatury
% figure(6);
% colormap(jet);
% contourf(X, R, transpose(energy.temp), 30);
% pbaspect([length radius 1]);
% 
% % Tempo reakcji SMR
% figure(7);
% colormap(jet);
% contourf(X,R,transpose(species_manager.R_st),30);
% pbaspect([length radius 1]);
% 
% % Tempo reakcji WGS
% figure(8);
% colormap(jet);
% contourf(X,R,transpose(species_manager.R_sh),30);
% pbaspect([length radius 1]);
% 
% % Średnie ułamki molowe wzdłuż reformera
% figure(9);
% species_list = values(species_manager.species);
% for i = 1:num_species
%     sp = species_list(i);
%     name = species_names(i);
%     molar_fractions = sp.Y ./ species_manager.Y_total .* species_manager.M_mix / sp.M;
%     avg_molar_fractions = sum(molar_fractions .* grid.face_area_x(1,:),2) / sum(grid.face_area_x(1,:),2);
%     plot(grid.cent_pos_x, avg_molar_fractions, 'DisplayName', name);
%     hold on;
% end
% 
% % Rezydua
% figure(10);
% semilogy(flow.residual_history_vx(1:flow.noi), 'DisplayName', 'vx residual');
% hold on;
% semilogy(flow.residual_history_vr(1:flow.noi), 'DisplayName', 'vr residual');
% semilogy(flow.residual_history_continuity(1:flow.noi), 'DisplayName', 'continuity residual');
% semilogy(energy.residual_history(1:energy.noi), 'DisplayName', 'T residual');
% for i = 1:num_species
%     name = species_names(i);
%     res = species_manager.residual_history(name);
%     semilogy(res(1:species_manager.noi), 'DisplayName', name + " residual");
% end

% Generacje entropii z każdego źródła:

ent_gen_chem_reac = entropy_generation_chemical_reaction(grid, energy, flow, species_manager);
ent_gen_heat_cond = entropy_generation_heat_transfer(grid, energy);
ent_gen_spc_diff = entropy_generation_species_diffusion(grid, energy, flow, species_manager);
ent_gen_flow_visc = entropy_generation_flow_viscosity(grid, energy, flow);

ent_gen_total = ent_gen_chem_reac + ent_gen_heat_cond + ent_gen_spc_diff + ent_gen_flow_visc;

entropy_generation = sum(sum(ent_gen_total .* grid.volume));
if crashed
    entropy_generation = Inf;
end

% Coefficients for interpolation
[coeff_x, coeff_r] = linear_interpolation_scheme(grid);

% Density
[coeff_x_temp, coeff_r_temp] = grid.domain_boundary.apply_boundary_condition_value(coeff_x, coeff_r);
[rho_x, ~] = evaluate_faces(coeff_x_temp, coeff_r_temp, flow.rho);

% CH4
[coeff_x_temp, coeff_r_temp] = ch4.species_bds.apply_boundary_condition_value(coeff_x, coeff_r);
[Y_ch4_x, ~] = evaluate_faces(coeff_x_temp, coeff_r_temp, ch4.Y);

% H2
[coeff_x_temp, coeff_r_temp] = h2.species_bds.apply_boundary_condition_value(coeff_x, coeff_r);
[Y_h2_x, ~] = evaluate_faces(coeff_x_temp, coeff_r_temp, h2.Y);

ch4_in = sum(Y_ch4_x(1,:) .* rho_x(1,:) .* grid.face_area_x(1,:) .* flow.vx_faces(1,:), 2);
ch4_out = sum(Y_ch4_x(end,:) .* rho_x(end,:) .* grid.face_area_x(end,:) .* flow.vx_faces(end,:), 2);
h2_out = sum(Y_h2_x(end,:) .* rho_x(end,:) .* grid.face_area_x(end,:) .* flow.vx_faces(end,:), 2);
end