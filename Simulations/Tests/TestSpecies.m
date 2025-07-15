%% Parametry symulacji
length = 0.3;
radius = 0.05;

porosity = 0.7;
permeability = 1e-7;
inertia_coefficient = 0.088;

solid_phase_thermal_conductivity = 20.0;
catalyst_density = 2.5e5;

temp_inlet = 800;
v_inlet = 1;
temp_wall = 800;
p_atm = 101325;
SC = 2;

a = 0.89;               % Order of reaction with respect to methane
b = 0.05;               % Order of reaction with respect to water
A_st = 1.354e-3;        % Arrhenius constant
E_a = 122500;           % Activation energy
delta_G = -28.6e3;      % Change of standard Gibbs free energy of water-gas-shift reaction
delta_H_st = 206e3;     % Enthalpy change accompanied with methane/steam reforming reaction
delta_H_sh = -41.2e3;   % Enthalpy change accompanied with water-gas-shift reaction

% Stała, której użycie wynika z użycia we wzorach na parametry jako
% argumentu temperatury podzielonej przez tysiąc
K = 1/1000;
% Stała w celu konwersji jednostek z uP na Pa*s
uP_to_Pas = 10^(-7);

% Parametry gazów
% Wodór
M_H2 = 2.016;
therm_cond_H2 = 0.01 * [-1.8954*K^6 11.972*K^5 -31.939*K^4 47.763*K^3 -47.190*K^2 62.892*K 1.5040];
visc_H2 = uP_to_Pas*[-9.9892*K^6 62.966*K^5 -167.51*K^4 249.41*K^3 -244.34*K^2 299.78*K 15.553];
heat_cap_H2 = [-6.4725*K^6 46.903*K^5 -136.15*K^4 199.29*K^3 -150.55*K^2 56.036*K 21.157];
diffusion_volume_H2 = 6.12;

% Tlenek węgla
M_CO = 28.010;
therm_cond_CO = 0.01 * [-2.3224*K^6 13.379*K^5 -30.818*K^4 36.018*K^3 -23.186*K^2 13.999*K -0.2815];
visc_CO = uP_to_Pas*[-32.298*K^6 208.42*K^5 -572.14*K^4 883.75*K^3 875.90*K^2 793.65*K -4.9137];
heat_cap_CO = [-7.6538*K^6 37.756*K^5 -66.346*K^4 41.974*K^3 5.2062*K^2 -8.1781*K 30.429];
diffusion_volume_CO = 18.0;

% Dwutlenek węgla
M_CO2 = 44.009;
therm_cond_CO2 = 0.01 * [18.698*K^6 -101.12*K^5 216.83*K^4 -233.29*K^3 129.65*K^2 -27.018*K 2.8888];
visc_CO2 = uP_to_Pas*[-0.4564*K^6 14.450*K^5 -85.929*K^4 244.22*K^3 -432.49*K^2 680.07*K -20.434];
heat_cap_CO2 = [-35.992*K^6 214.58*K^5 -519.9*K^4 657.88*K^3 -471.33*K^2 204.60*K 4.3669];
diffusion_volume_CO2 = 26.7;

% Metan
M_CH4 = 16.043;
therm_cond_CH4 = 0.01 * [3.2774*K^6 -17.283*K^5 38.251*K^4 -47.440*K^3 37.413*K^2 1.8732*K 0.4796];
visc_CH4 = uP_to_Pas*[-22.920*K^6 140.48*K^5 -367.06*K^4 548.11*K^3 -543.82*K^2 529.37*K -9.9989];
heat_cap_CH4 = [61.321*K^6 -358.75*K^5 856.93*K^4 -1068.7*K^3 712.55*K^2 -178.59*K 47.964];
diffusion_volume_CH4 = 25.14;

% Para wodna
M_H2O = 18.015;
therm_cond_H2O = 0.01 * [4.1531*K^6 -18.974*K^5 35.993*K^4 -41.390*K^3 35.922*K^2 -7.9139*K 2.0103];
visc_H2O = uP_to_Pas*[19.591*K^6 -126.96*K^5 348.12*K^4 -522.38*K^3 419.50*K^2 244.93*K -6.7541];
heat_cap_H2O = [14.015*K^6 -79.409*K^5 181.54*K^4 -217.08*K^3 146.01*K^2 -41.205*K 37.373];
diffusion_volume_H2O = 13.1;

% Parametry solvera
inner_iters = 5;
rel_fact_sp = 0.75;
rel_fact_R_st = 0.3;
rel_fact_R_sh = 0.1;
tol_sp = 1e-7;
sp_iters = 25;

%% Definicja siatki
sz = [60 10];

offset = 1e-30;
dx = length * ones([sz(1) 1]) / sz(1);
dr = radius * ones([1 sz(2)]) / sz(2);

grid = Grid2D(dx, dr, offset);

%% Definicja komponentów

flow = FlowComponent(grid, v_inlet*ones(sz), zeros(sz), p_atm * ones(sz), 1, @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), 1e-6, 1e-6, 0.55, 0.2, 1, 1, 30);
flow.vx_faces = v_inlet*ones(sz + [1 0]);
flow.vr_faces = zeros(sz + [0 1]);
energy = EnergyComponent(grid, flow, 1200*ones(sz), @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), 1e-9, 1, 5, 30);

species_manager = SpeciesManagerComponent(grid, flow, energy, @(x) 0, a, b, A_st, E_a, delta_G, catalyst_density, rel_fact_R_st, rel_fact_R_sh, inner_iters);
h2 = Species("H2", grid, flow, species_manager, zeros(sz), heat_cap_H2, visc_H2, therm_cond_H2, M_H2, diffusion_volume_H2, tol_sp, rel_fact_sp, sp_iters);
co = Species("CO", grid, flow, species_manager, zeros(sz), heat_cap_CO, visc_CO, therm_cond_CO, M_CO, diffusion_volume_CO, tol_sp, rel_fact_sp, sp_iters);
co2 = Species("CO2", grid, flow, species_manager, zeros(sz), heat_cap_CO2, visc_CO2, therm_cond_CO2, M_CO2, diffusion_volume_CO2, tol_sp, rel_fact_sp, sp_iters);
ch4 = Species("CH4", grid, flow, species_manager, M_CH4/(M_CH4+SC*M_H2O) * ones(sz), heat_cap_CH4, visc_CH4, therm_cond_CH4, M_CH4, diffusion_volume_CH4, tol_sp, rel_fact_sp, sp_iters);
h2o = Species("H2O", grid, flow, species_manager, SC*M_H2O/(M_CH4+SC*M_H2O) * ones(sz), heat_cap_H2O, visc_H2O, therm_cond_H2O, M_H2O, diffusion_volume_H2O, tol_sp, rel_fact_sp, sp_iters);

properties_manager = ThermophysicalProperties(grid, flow, energy, species_manager, 20);
species_manager.D_f = @properties_manager.D_function;

species_manager.update_properties();
flow.rho_function = @properties_manager.rho_function;
flow.update_properties();
%% Warunki brzegowe
inlet_ch4 = FixedValueBoundary("Inlet CH4", grid, M_CH4/(M_CH4+SC*M_H2O));
inlet_h2o = FixedValueBoundary("Inlet H2O", grid, SC*M_H2O/(M_CH4+SC*M_H2O));
inlet_rest = FixedValueBoundary("Inlet rest", grid, 0);
wall_and_outlet_all = FixedNormalDerivativeBoundary("Rest", grid, 0);

domain_boundary = FixedNormalDerivativeBoundary("Domain", grid, 0);


for j = 1:sz(2)
    % Inlet
    inlet_ch4.add_boundary_face("x", [1 j], 1);
    inlet_h2o.add_boundary_face("x", [1 j], 1);
    inlet_rest.add_boundary_face("x", [1 j], 1);
    domain_boundary.add_boundary_face("x", [1 j], 1);

    % Outlet
    wall_and_outlet_all.add_boundary_face("x", [sz(1)+1 j], -1);
    domain_boundary.add_boundary_face("x", [sz(1)+1 j], -1)
end

for i = 1:sz(1)
    % Axis
    wall_and_outlet_all.add_boundary_face("r", [i 1], 1);
    domain_boundary.add_boundary_face("r", [i 1], 1);

    % Wall
    wall_and_outlet_all.add_boundary_face("r", [i sz(2)+1], -1);
    domain_boundary.add_boundary_face("r", [i sz(2)+1], -1);
end

grid.domain_boundary.add_boundary(domain_boundary);
h2.species_bds.add_boundary(inlet_rest);
h2.species_bds.add_boundary(wall_and_outlet_all);
co.species_bds.add_boundary(inlet_rest);
co.species_bds.add_boundary(wall_and_outlet_all);
co2.species_bds.add_boundary(inlet_rest);
co2.species_bds.add_boundary(wall_and_outlet_all);
ch4.species_bds.add_boundary(inlet_ch4);
ch4.species_bds.add_boundary(wall_and_outlet_all);
h2o.species_bds.add_boundary(inlet_h2o);
h2o.species_bds.add_boundary(wall_and_outlet_all);


names = ["H2", "CO", "CO2", "CH4", "H2O"];
a = false;
while ~a
    a = species_manager.iterate();
    if mod(species_manager.noi,5)==0
        fprintf('Current iteration: %d \n', species_manager.noi);
        for i = 1:5
            res_hist = species_manager.residual_history(names(i));
            fprintf('The %s residual is equal to: %d\n', names(i), res_hist(species_manager.noi));
        end
    end
end

species_list = values(species_manager.species);
species_names = keys(species_manager.species);
num_species = numel(species_list);

% Grid
[X,R]=meshgrid(grid.cent_pos_x, grid.cent_pos_r);

% Values have to be transposed - probably because MATLAB uses column major order
% Steam-methane reforming reactor rate
figure(1);
colormap(jet);
contourf(X,R,transpose(species_manager.R_st),30);

% Water-gas-shift reaction rate
figure(2);
colormap(jet);
contourf(X,R,transpose(species_manager.R_sh),30);

% Average molar fractions along the reformer
figure(3);
for i = 1:num_species
    sp = species_list(i);
    name = names(i);
    molar_fractions = sp.Y ./ species_manager.Y_total .* species_manager.M_mix / sp.M;
    avg_molar_fractions = sum(molar_fractions .* grid.face_area_x(1,:),2) / sum(grid.face_area_x(1,:),2);
    plot(grid.cent_pos_x, avg_molar_fractions, 'DisplayName', name);
    hold on;
end

% Residuals
figure(4);
for i = 1:num_species
    name = names(i);
    res = species_manager.residual_history(name);
    semilogy(res(1:species_manager.noi),'DisplayName',name + " residual");
    hold on;
end
