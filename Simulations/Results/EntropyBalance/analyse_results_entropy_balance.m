path = "D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\GridSensitivityAnalysis\saved_results\";
image_path = "C:\Users\zgolp\OneDrive\Pulpit\AGH\INŻYNIERKA\ZDJĘCIA\BILANS\";
reference_porosity = 0.65;
reference_catalyst_density = 2.5e5;
reference_wall_temp = 900;
reference_inlet_v = 0.25;
reference_SC = 2;

grid_sizes = (2:12) / 2;

ent_gen_tot = zeros(size(grid_sizes));
ent_gen_tot_balance = zeros(size(grid_sizes));
ent_gen_rel_diff = zeros(size(grid_sizes));
ent_flux_inlet = zeros(size(grid_sizes));
ent_flux_outlet = zeros(size(grid_sizes));
ent_flux_wall = zeros(size(grid_sizes));
ent_flux_inlet_conv = zeros(size(grid_sizes));
ent_flux_inlet_heat = zeros(size(grid_sizes));
ent_flux_inlet_diff = zeros(size(grid_sizes));

for i = 1:numel(grid_sizes)
    S = load(get_file_path_grid_analysis(grid_sizes(i), reference_porosity, reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, path));
    ent_gen_tot(i) = S.ent_gen_total;

    [ent_gen_rel_diff(i), ent_gen_tot_balance(i), ent_flux_inlet(i), ent_flux_outlet(i), ent_flux_wall(i), ent_flux_inlet_conv(i), ent_flux_inlet_heat(i), ent_flux_inlet_diff(i)] = entropy_balance(S.grid, S.energy, S.flow, S.species_manager, S.ent_gen_total);
end


figure; 
set(gcf,'position', [100 100 600 600]);
pbaspect([6 4 1]);
hold on;
box on;
grid_sizes = 2 * grid_sizes;
RGB = orderedcolors("gem");
RGB(8,:) = [0.2240 0.5540 0.5040];
colororder(RGB);

set(gca,'fontsize', 11.5);
ylabel('$$\Delta_\mathrm{rel} \dot{S}_\mathrm{gen} \ [-]$$', 'Interpreter', 'latex', 'FontSize', 18);
xlabel("$$n \ [-]$$", 'Interpreter', 'latex', 'FontSize', 18),
scatter(grid_sizes, ent_gen_rel_diff, 400, ".");
xlim([1.7 12.3]);
ylim([-0.14, 0]);
xticks(2:12);
ax = gca;
yticks(-0.14:0.02:0);
ax.YTickLabel = {'-14%', '-12%', '-10%', '-8%', '-6%', '-4%', '-2%', '0%'};

exportgraphics(gcf, image_path + 'entropy_balance_relative_difference.pdf');