offset_R = 1e-30;
delta_R = 0.05;
L = 0.3;

sz = [90 15];

dx = L * ones(sz(1),1) / sz(1);
dr = delta_R * ones(1,sz(2)) / sz(2);

grid = Grid2D(dx, dr, offset_R);

flow = FlowComponent(grid, zeros(sz), zeros(sz), zeros(sz), 1, @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), 1e-7, 1e-7, 0.7, 0.25, 1, 1, 70);

% Boundary conditions
% West
for j = 1:sz(2)
    flow.vx_bds.set_boundary_condition("x", [1 j], 1, 1, 0, 1);
    flow.vr_bds.set_boundary_condition("x", [1 j], 1, 1, 0, 0);
    flow.p_bds.set_boundary_condition("x", [1 j], 1, 0, 1, 0);

    grid.domain_boundary.set_boundary_condition("x", [1 j], 1, 0, 1, 0);
end
% East
for j = 1:sz(2)
    flow.vx_bds.set_boundary_condition("x", [sz(1)+1 j], -1, 0, 1, 0);
    flow.vr_bds.set_boundary_condition("x", [sz(1)+1 j], -1, 0, 1, 0);
    flow.p_bds.set_boundary_condition("x", [sz(1)+1 j], -1, 1, 0, 0);

    grid.domain_boundary.set_boundary_condition("x", [sz(1)+1 j], -1, 0, 1, 0);
end
% South
for i = 1:sz(1)
    flow.vx_bds.set_boundary_condition("r", [i 1], 1, 0, 1, 0);
    flow.vr_bds.set_boundary_condition("r", [i 1], 1, 1, 0, 0);
    flow.p_bds.set_boundary_condition("r", [i 1], 1, 0, 1, 0);

    grid.domain_boundary.set_boundary_condition("r", [i 1], 1, 0, 1, 0);
end
% North
for i = 1:sz(1)
    flow.vx_bds.set_boundary_condition("r", [i sz(2)+1], -1, 1, 0, 0);
    flow.vr_bds.set_boundary_condition("r", [i sz(2)+1], -1, 1, 0, 0);
    flow.p_bds.set_boundary_condition("r", [i sz(2)+1], -1, 0, 1, 0);

    grid.domain_boundary.set_boundary_condition("r", [i sz(2)+1], -1, 0, 1, 0);
end

flow.convert_pressure_bd_conditions();

flow.update_properties();
[coeff_vx, coeff_vr] = flow.get_coefficients_v();
flow.update_face_velocities(flow.vx, flow.vr, coeff_vx(:,:,1), coeff_vr(:,:,1));

a = false;
while ~a
    a = flow.iterate();
    if mod(flow.noi,10)==0
        fprintf('Current iteration: %d \n', flow.noi);
        fprintf('The vx residual is equal to: %d \n', flow.residual_history_vx(flow.noi));
        fprintf('The vr residual is equal to: %d \n', flow.residual_history_vr(flow.noi));
        fprintf('The p residual is equal to: %d \n', flow.residual_history_p(flow.noi));
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
semilogy(flow.residual_history_p(1:flow.noi),'DisplayName','p residual');