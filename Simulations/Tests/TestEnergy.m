offset_R = 1;
delta_R = 1;
L = 1;

sz = [50 50];

dx = L * ones(sz(1),1) / sz(1);
dr = delta_R * ones(1,sz(2)) / sz(2);

grid = Grid2D(dx, dr, offset_R);

rho = 1000; % Density [kg/m^3]
visc = 1; % Dynamic viscosity [Pa*s]

vx = 0.01; % Velocity in the x direction [m/s]
vx_field = vx * ones(sz); % Velocity field in the x direction [m/s]
vr_field = zeros(sz); % Velocity field in the r direction [m/s]
p = zeros(sz); % Pressure [Pa]
porosity = ones(sz); % Porosity [-]
rho_function = @(x) rho * ones(sz); % Density field [kg/m^3]
visc_function = @(x) visc * ones(sz); % Dynamic viscosity field [Pa*s]
source_terms_flow = @(x) zeros(sz); % Source terms of the momentum equation - set to zero

T_init = zeros(sz); % Initial temperature field [K]
k = @(x) ones(sz); % Thermal conductivity [W/(m*K)]
Cp = @(x) ones(sz); % Heat capacity [J/(kg*K)]
q = @(x) ones(sz); % Heat source [W/m^3]
q_t = @(x) zeros(sz); % Linear heat source [W/(m^3*K)]
rel_fact_T = 0.9; 
tol_T = 1e-7;
inner_iters = 1;
solver_iters = 10;


flow = FlowComponent(grid, vx_field, vr_field, p, porosity, rho_function, visc_function, source_terms_flow, source_terms_flow, source_terms_flow, source_terms_flow, 1, 1, 1, 1, 1, 1, 1); % Relaxation factors are set to 1 because the component is never updated
energy = EnergyComponent(grid, flow, T_init, k, Cp, q, q_t, tol_T, rel_fact_T, inner_iters, solver_iters);

% Boundary conditions
T_inlet_outlet = 100; % Temperature at the inlet and the outlet [K]
T_outer_wall = 200; % Outer wall temperature [K]
T_inner_wall = 400; % Inner wall temperature [K]
inlet_vx = FixedValueBoundary("Inlet vx", grid, vx);
inlet_vr = FixedValueBoundary("Inlet vr", grid, 0);
inlet_p = FixedNormalDerivativeBoundary("Inlet p", grid, 0);
inlet_temp = FixedValueBoundary("Inlet T", grid, T_inlet_outlet);

outlet_vx = FixedValueBoundary("Outlet vx", grid, vx);
outlet_vr = FixedValueBoundary("Outlet vr", grid, 0);
outlet_p = FixedNormalDerivativeBoundary("Outlet p", grid, 0);
outlet_temp = FixedValueBoundary("Outlet T", grid, T_inlet_outlet);

outer_wall_vx = FixedValueBoundary("Outer wall vx", grid, 0);
outer_wall_vr = FixedValueBoundary("Outer wall vr", grid, 0);
outer_wall_p = FixedNormalDerivativeBoundary("Outer wall p", grid, 0);
outer_wall_temp = FixedValueBoundary("Outer wall T", grid, T_outer_wall);

inner_wall_vx = FixedValueBoundary("Inner wall vx", grid, 0);
inner_wall_vr = FixedValueBoundary("Inner wall vr", grid, 0);
inner_wall_p = FixedNormalDerivativeBoundary("Inner wall p", grid, 0);
inner_wall_temp = FixedValueBoundary("Inner wall T", grid, T_inner_wall);

domain_boundary = FixedNormalDerivativeBoundary("Domain", grid, 0);


for j = 1:sz(2)
    % Inlet
    inlet_vx.add_boundary_face("x", [1 j], 1);
    inlet_vr.add_boundary_face("x", [1 j], 1);
    inlet_p.add_boundary_face("x", [1 j], 1);
    inlet_temp.add_boundary_face("x", [1 j], 1);
    domain_boundary.add_boundary_face("x", [1 j], 1);

    % Outlet
    outlet_vx.add_boundary_face("x", [sz(1)+1 j], -1);
    outlet_vr.add_boundary_face("x", [sz(1)+1 j], -1);
    outlet_p.add_boundary_face("x", [sz(1)+1 j], -1);
    outlet_temp.add_boundary_face("x", [sz(1)+1 j], -1);
    domain_boundary.add_boundary_face("x", [sz(1)+1 j], -1)
end

% Axis
for i = 1:sz(1)
    % Inner wall
    inner_wall_vx.add_boundary_face("r", [i 1], 1);
    inner_wall_vr.add_boundary_face("r", [i 1], 1);
    inner_wall_p.add_boundary_face("r", [i 1], 1);
    inner_wall_temp.add_boundary_face("r", [i 1], 1);
    domain_boundary.add_boundary_face("r", [i 1], 1);

    % Outer wall
    outer_wall_vx.add_boundary_face("r", [i sz(2)+1], -1);
    outer_wall_vr.add_boundary_face("r", [i sz(2)+1], -1);
    outer_wall_p.add_boundary_face("r", [i sz(2)+1], -1);
    outer_wall_temp.add_boundary_face("r", [i sz(2)+1], -1);
    domain_boundary.add_boundary_face("r", [i sz(2)+1], -1);
end

grid.domain_boundary.add_boundary(domain_boundary);

flow.vx_bds.add_boundary(inlet_vx);
flow.vx_bds.add_boundary(outlet_vx);
flow.vx_bds.add_boundary(inner_wall_vx);
flow.vx_bds.add_boundary(outer_wall_vx);

flow.vr_bds.add_boundary(inlet_vr);
flow.vr_bds.add_boundary(outlet_vr);
flow.vr_bds.add_boundary(inner_wall_vr);
flow.vr_bds.add_boundary(outer_wall_vr);

flow.p_bds.add_boundary(inlet_p);
flow.p_bds.add_boundary(outlet_p);
flow.p_bds.add_boundary(inner_wall_p);
flow.p_bds.add_boundary(outer_wall_p);

energy.temp_bds.add_boundary(inlet_temp);
energy.temp_bds.add_boundary(outlet_temp);
energy.temp_bds.add_boundary(inner_wall_temp);
energy.temp_bds.add_boundary(outer_wall_temp);

flow.convert_pressure_bd_conditions();

flow.update_properties();
[coeff_vx, coeff_vr] = flow.get_coefficients_v();
flow.update_face_velocities(flow.vx, flow.vr, coeff_vx(:,:,1), coeff_vr(:,:,1));
clear coeff_vx coeff_vr

energy.q_function = @(x) 10*ones(sz);
energy.update_properties();

a = false;
while ~a
    a = energy.iterate();
    if mod(energy.noi,10)==0
        fprintf('Current iteration: %d \n', energy.noi);
        fprintf('The T residual is equal to: %d \n', energy.residual_history(energy.noi));
    end
end


% Grid
[X,R]=meshgrid(grid.cent_pos_x, grid.cent_pos_r);

% Values have to be transposed
% Temperature graph
figure(1);
colormap(jet);
contourf(X,R,transpose(energy.temp),30);
pbaspect([L delta_R, 1]);

% Residuals
figure(2);
semilogy(energy.residual_history(1:energy.noi),'DisplayName','T residual');