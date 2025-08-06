offset_R = 1e-30;
delta_R = 0.05;
L = 0.3;

sz = [60 10];

dx = L * ones(sz(1),1) / sz(1);
dr = delta_R * ones(1,sz(2)) / sz(2);

grid = Grid2D(dx, dr, offset_R);

flow = FlowComponent(grid, zeros(sz), zeros(sz), zeros(sz), 1, @(x) 1 * ones(sz), @(x) 1 * ones(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), 1e-7, 1e-7, 0.7, 0.3, 1, 1, 70);

% Boundary conditions
inlet_vx = FixedValueBoundary("Inlet vx", grid, 0.01);
inlet_vr = FixedValueBoundary("Inlet vr", grid, 0);
inlet_p = FixedNormalDerivativeBoundary("Inlet p", grid, 0);
% inlet_p = FixedNormalDerivativeBoundary("Inlet p", grid, 0.01*(-1 / 1e-7 - 1 * 0.088 / sqrt(1e-7) .* sqrt(0.01 .^ 2))); % Use with the additional source terms

outlet_vx = FixedNormalDerivativeBoundary("Outlet vx", grid, 0);
outlet_vr = FixedNormalDerivativeBoundary("Outlet vr", grid, 0);
outlet_p = FixedValueBoundary("Outlet p", grid, 0);

wall_vx = FixedValueBoundary("Wall vx", grid, 0);
wall_vr = FixedValueBoundary("Wall vr", grid, 0);
wall_p = FixedNormalDerivativeBoundary("Wall p", grid, 0);

axis_vx = FixedNormalDerivativeBoundary("Axis vx", grid, 0);
axis_vr = FixedValueBoundary("Axis vr", grid, 0);
axis_p = FixedNormalDerivativeBoundary("Axis p", grid, 0);

domain_boundary = FixedNormalDerivativeBoundary("Domain", grid, 0);


for j = 1:sz(2)
    % Inlet
    inlet_vx.add_boundary_face("x", [1 j], 1);
    inlet_vr.add_boundary_face("x", [1 j], 1);
    inlet_p.add_boundary_face("x", [1 j], 1);
    domain_boundary.add_boundary_face("x", [1 j], 1);

    % Outlet
    outlet_vx.add_boundary_face("x", [sz(1)+1 j], -1);
    outlet_vr.add_boundary_face("x", [sz(1)+1 j], -1);
    outlet_p.add_boundary_face("x", [sz(1)+1 j], -1);
    domain_boundary.add_boundary_face("x", [sz(1)+1 j], -1)
end

% Axis
for i = 1:sz(1)
    % Axis
    axis_vx.add_boundary_face("r", [i 1], 1);
    axis_vr.add_boundary_face("r", [i 1], 1);
    axis_p.add_boundary_face("r", [i 1], 1);
    domain_boundary.add_boundary_face("r", [i 1], 1);

    % Wall
    wall_vx.add_boundary_face("r", [i sz(2)+1], -1);
    wall_vr.add_boundary_face("r", [i sz(2)+1], -1);
    wall_p.add_boundary_face("r", [i sz(2)+1], -1);
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

% flow.src_linx_function = @() lin_src_vx(flow, 1e-7, 0.088);
% flow.src_linr_function = @() lin_src_vr(flow, 1e-7, 0.088);

flow.update_properties();
[coeff_vx, coeff_vr] = flow.get_coefficients_v();
flow.update_face_velocities(flow.vx, flow.vr, coeff_vx(:,:,1), coeff_vr(:,:,1));

a = false;
while ~a
    a = flow.iterate();
    if mod(flow.noi,1)==0
        fprintf('Current iteration: %d \n', flow.noi);
        fprintf('The vx residual is equal to: %d \n', flow.residual_history_vx(flow.noi));
        fprintf('The vr residual is equal to: %d \n', flow.residual_history_vr(flow.noi));
        fprintf('The continuity residual is equal to: %d \n', flow.residual_history_continuity(flow.noi));
    end
end

% Grid
[X, R]=meshgrid(grid.cent_pos_x, grid.cent_pos_r);


% Values have to be transposed - probably because MATLAB uses column major order
 
% Velocity magnitude graph
figure(1);
colormap(jet);
contourf(X,R,transpose(sqrt(flow.vx .^ 2 + flow.vr .^ 2)),30);

% Velocity in the X direction graph
figure(2);
colormap(jet);
contourf(X,R,transpose(flow.vx),30);

% Velocity in the R direction graph
figure(3);
colormap(jet);
contourf(X,R,transpose(flow.vr),30);

% Pressure graph
figure(4);
colormap(jet);
contourf(X,R,transpose(flow.p),30);

% Streamlines
figure(5);
streamslice(X,R,transpose(flow.vx),transpose(flow.vr),2);

% Residuals
figure(6);
semilogy(flow.residual_history_vx(1:flow.noi),'DisplayName','vx residual');
hold on;
semilogy(flow.residual_history_vr(1:flow.noi),'DisplayName','vr residual');
semilogy(flow.residual_history_continuity(1:flow.noi),'DisplayName','continuity residual');