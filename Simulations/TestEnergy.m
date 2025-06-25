offset_R = 0;
delta_R = 1;
L = 1;

sz = [50 50];

dx = L * ones(sz(1),1) / sz(1);
dr = delta_R * ones(1,sz(2)) / sz(2);

grid = Grid2D(dx, dr, offset_R);

flow = FlowComponent(grid, 1*ones(sz), zeros(sz), zeros(sz), 1, @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), 1e-6, 1e-6, 0.55, 0.2, 1, 1, 30);
energy = EnergyComponent(grid, flow, zeros(sz), @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), 1e-9, 1, 5, 30);

% Boundary conditions
% West
for j = 1:sz(2)
    flow.vx_bds.set_boundary_condition("x", [1 j], 1, 1, 0, 1);
    flow.vr_bds.set_boundary_condition("x", [1 j], 1, 1, 0, 0);
    flow.p_bds.set_boundary_condition("x", [1 j], 1, 1, 0, 0);

    energy.temp_bds.set_boundary_condition("x", [1 j], 1, 1, 0, 100);
end
% East
for j = 1:sz(2)
    flow.vx_bds.set_boundary_condition("x", [sz(1)+1 j], -1, 1, 0, 1);
    flow.vr_bds.set_boundary_condition("x", [sz(1)+1 j], -1, 1, 0, 0);
    flow.p_bds.set_boundary_condition("x", [sz(1)+1 j], -1, 1, 0, 0);

    energy.temp_bds.set_boundary_condition("x", [sz(1)+1 j], -1, 0, 1, 0);
end
% South
for i = 1:sz(1)
    flow.vx_bds.set_boundary_condition("r", [i 1], 1, 1, 0, 1);
    flow.vr_bds.set_boundary_condition("r", [i 1], 1, 1, 0, 0);
    flow.p_bds.set_boundary_condition("r", [i 1], 1, 1, 0, 0);

    energy.temp_bds.set_boundary_condition("r", [i 1], 1, 0, 1, 0);
end
% North
for i = 1:sz(1)
    flow.vx_bds.set_boundary_condition("r", [i sz(2)+1], -1, 1, 0, 1);
    flow.vr_bds.set_boundary_condition("r", [i sz(2)+1], -1, 1, 0, 0);
    flow.p_bds.set_boundary_condition("r", [i sz(2)+1], -1, 1, 0, 0);

    energy.temp_bds.set_boundary_condition("r", [i sz(2)+1], -1, 1, 0, 300);
end

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