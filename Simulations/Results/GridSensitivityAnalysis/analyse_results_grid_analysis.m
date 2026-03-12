path = "D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\GridSensitivityAnalysis\saved_results\";
image_path = "C:\Users\zgolp\OneDrive\Pulpit\AGH\INŻYNIERKA\ZDJĘCIA\SIATKA\";
reference_porosity = 0.65;
reference_catalyst_density = 2.5e5;
reference_wall_temp = 900;
reference_inlet_v = 0.25;
reference_SC = 2;

grid_sizes = (2:12) / 2;

ent_gen_tot = zeros(size(grid_sizes));
ent_gen_CR = zeros(size(grid_sizes));
ent_gen_HT = zeros(size(grid_sizes));
ent_gen_FV = zeros(size(grid_sizes));
ent_gen_SD = zeros(size(grid_sizes));

h2_prod = zeros(size(grid_sizes));
mcr = zeros(size(grid_sizes));
spec_ent_gen = zeros(size(grid_sizes));

ent_gen_rel_diff = zeros(size(grid_sizes));
ent_gen_tot_balance = zeros(size(grid_sizes));
ent_flux_inlet = zeros(size(grid_sizes));
ent_flux_outlet = zeros(size(grid_sizes));
ent_flux_wall = zeros(size(grid_sizes));
ent_flux_inlet_conv = zeros(size(grid_sizes));
ent_flux_inlet_heat = zeros(size(grid_sizes));
ent_flux_inlet_diff = zeros(size(grid_sizes));

for i = 1:numel(grid_sizes)
    S = load(get_file_path_grid_analysis(grid_sizes(i), reference_porosity, reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, path));
    ent_gen_tot(i) = S.ent_gen_total;
    ent_gen_CR(i) = S.ent_gen_whole_reformer("Chemical reactions");
    ent_gen_HT(i) = S.ent_gen_whole_reformer("Heat transfer");
    ent_gen_FV(i) = S.ent_gen_whole_reformer("Flow viscosity");
    ent_gen_SD(i) = S.ent_gen_whole_reformer("Species diffusion");

    h2_prod(i) = S.m_h2_out;
    mcr(i) = S.methane_conversion_rate;
    spec_ent_gen(i) = S.entropy_per_hydrogen;

    [ent_gen_rel_diff(i), ent_gen_tot_balance(i), ent_flux_inlet(i), ent_flux_outlet(i), ent_flux_wall(i), ent_flux_inlet_conv(i), ent_flux_inlet_heat(i), ent_flux_inlet_diff(i)] = entropy_balance(S.grid, S.energy, S.flow, S.species_manager, S.ent_gen_total);
end

figure; 
set(gcf,'position',[100 100 600 400]);
hold on;
box on;
grid_sizes = 2 * grid_sizes;
RGB = orderedcolors("gem");
RGB(8,:) = [0.2240 0.5540 0.5040];
colororder(RGB);

s1 = scatter(grid_sizes, ent_gen_tot / ent_gen_tot(end), 400, ".");
s2 = scatter(grid_sizes, ent_gen_CR / ent_gen_CR(end), 400, ".");
s3 = scatter(grid_sizes, ent_gen_HT / ent_gen_HT(end), 400, ".");
s4 = scatter(grid_sizes, ent_gen_FV / ent_gen_FV(end), 400, ".");
s5 = scatter(grid_sizes, ent_gen_SD / ent_gen_SD(end), 400, ".");
s6 = scatter(grid_sizes, h2_prod / h2_prod(end), 400, ".");
s7 = scatter(grid_sizes, mcr / mcr(end), 400, ".");
s8 = scatter(grid_sizes, spec_ent_gen / spec_ent_gen(end), 400, ".");

% Create legend with dummy markers to control exact size
h = zeros(1,8);
h(1) = plot(nan, nan, '.', 'MarkerSize', 18, 'Color', s1.CData(1,:));
h(2) = plot(nan, nan, '.', 'MarkerSize', 18, 'Color', s2.CData(1,:));
h(3) = plot(nan, nan, '.', 'MarkerSize', 18, 'Color', s3.CData(1,:));
h(4) = plot(nan, nan, '.', 'MarkerSize', 18, 'Color', s4.CData(1,:));
h(5) = plot(nan, nan, '.', 'MarkerSize', 18, 'Color', s5.CData(1,:));
h(6) = plot(nan, nan, '.', 'MarkerSize', 18, 'Color', s6.CData(1,:));
h(7) = plot(nan, nan, '.', 'MarkerSize', 18, 'Color', s7.CData(1,:));
h(8) = plot(nan, nan, '.', 'MarkerSize', 18, 'Color', s8.CData(1,:));

lgd = legend(h, ...
    '$$\dot{S}_\mathrm{gen}$$', '$$\dot{S}_\mathrm{gen,CR}$$', ...
    '$$\dot{S}_\mathrm{gen,HT}$$', '$$\dot{S}_\mathrm{gen,FV}$$', ...
    '$$\dot{S}_\mathrm{gen,SD}$$', '$$\dot{m}_\mathrm{H_2,out}$$', ...
    '$$\mathit{\Phi}_\mathrm{CH_4}$$', '$$s_\mathrm{gen,H_2}$$', ...
    "Location", "southeast", 'FontSize', 13);

set(lgd, 'Interpreter', 'latex');

set(lgd, 'Interpreter', 'latex');
set(gca,'fontsize', 11.5);
set(ylabel("$$y_n / y_{12} \ [-]$$", 'FontSize', 18), 'Interpreter', 'latex');
set(xlabel("$$n \ [-]$$", 'FontSize', 18), 'Interpreter', 'latex');
xticks(2:12);
axis([1.7 12.3 0.7 1.05]);


exportgraphics(gcf, image_path + 'grid_sensitivity_analysis.pdf');


% Pole generacji ze względu na dyfuzję
S = load(get_file_path_grid_analysis(grid_sizes(11) / 2, reference_porosity, reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, path));
radius = 0.05;
length = 0.3;
[X, R] = meshgrid(S.grid.cent_pos_x, S.grid.cent_pos_r);
X = (X - X(1)) / (X(end) - X(1)) * length / radius;
R = (R - R(1)) / (R(end) - R(1));


figure;
contourf(X, R, transpose(S.ent_gen("Species diffusion")), 10);
colormap(jet);

set(gcf, 'position', [100 100 1200 300]);
xlim([0 6]);
ylim([0 1]);
cb = colorbar('FontSize', 16);

ax = gca;
ax.YAxis.FontSize = 16;
ax.XAxis.FontSize = 16;

set(ylabel("$$r/R \ [-]$$", 'FontSize', 20), 'Interpreter', 'latex');
set(xlabel("$$x/R \ [-]$$", 'FontSize', 20), 'Interpreter', 'latex');
set(title("$$\sigma_\mathrm{SD} \ [\mathrm{W \, m^{-3} \, K^{-1}}]$$", 'FontSize', 23), 'Interpreter', 'latex');
axis equal;

exportgraphics(gcf, image_path + 'grid_sensitivity_analysis_sd_contour.pdf');