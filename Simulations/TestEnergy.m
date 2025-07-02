offset_R = 0;
delta_R = 1;
L = 1;

sz = [50 50];

dx = L * ones(sz(1),1) / sz(1);
dr = delta_R * ones(1,sz(2)) / sz(2);

grid = Grid2D(dx, dr, offset_R);

flow = FlowComponent(grid, 5*ones(sz), zeros(sz), zeros(sz), 1, @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), 1e-6, 1e-6, 0.55, 0.2, 1, 1, 30);
energy = EnergyComponent(grid, flow, zeros(sz), @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), 1e-9, 1, 5, 30);

% Boundary conditions
inlet_vx = FixedValueBoundary("Inlet vx", grid, 5);
inlet_vr = FixedValueBoundary("Inlet vr", grid, 0);
inlet_p = FixedNormalDerivativeBoundary("Inlet p", grid, 0);
inlet_temp = FixedValueBoundary("Inlet T", grid, 100);

outlet_vx = FixedValueBoundary("Outlet vx", grid, 5);
outlet_vr = FixedValueBoundary("Outlet vr", grid, 0);
outlet_p = FixedNormalDerivativeBoundary("Outlet p", grid, 0);
outlet_temp = FixedValueBoundary("Outlet T", grid, 100);

wall_vx = FixedValueBoundary("Wall vx", grid, 0);
wall_vr = FixedValueBoundary("Wall vr", grid, 0);
wall_p = FixedNormalDerivativeBoundary("Wall p", grid, 0);
wall_temp = FixedValueBoundary("Wall T", grid, 200);

axis_vx = FixedValueBoundary("Axis vx", grid, 0);
axis_vr = FixedValueBoundary("Axis vr", grid, 0);
axis_p = FixedNormalDerivativeBoundary("Axis p", grid, 0);
axis_temp = FixedValueBoundary("Axis T", grid, 400);

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
    % Axis
    axis_vx.add_boundary_face("r", [i 1], 1);
    axis_vr.add_boundary_face("r", [i 1], 1);
    axis_p.add_boundary_face("r", [i 1], 1);
    axis_temp.add_boundary_face("r", [i 1], 1);
    domain_boundary.add_boundary_face("r", [i 1], 1);

    % Wall
    wall_vx.add_boundary_face("r", [i sz(2)+1], -1);
    wall_vr.add_boundary_face("r", [i sz(2)+1], -1);
    wall_p.add_boundary_face("r", [i sz(2)+1], -1);
    wall_temp.add_boundary_face("r", [i sz(2)+1], -1);
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

energy.temp_bds.add_boundary(inlet_temp);
energy.temp_bds.add_boundary(outlet_temp);
energy.temp_bds.add_boundary(axis_temp);
energy.temp_bds.add_boundary(wall_temp);

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

% Values have to be transposed - probably because MATLAB uses column major order
% Temperature graph
figure(1);
colormap(jet);
contourf(X,R,transpose(energy.temp),30);

% Residuals
figure(2);
semilogy(energy.residual_history(1:energy.noi),'DisplayName','T residual');