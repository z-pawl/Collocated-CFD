case_data = load('D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\OneChanged\saved_results\reference_case.mat');
radius = 0.05;
length = 0.3;
[X, R] = meshgrid(case_data.grid.cent_pos_x, case_data.grid.cent_pos_r);
X = (X - X(1)) / (X(end) - X(1)) * length / radius;
R = (R - R(1)) / (R(end) - R(1));

figure;

case_data = load('D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\OneChanged\saved_results\reference_case.mat');
Y1 = transpose((case_data.ent_gen("Chemical reactions") + ...
    case_data.ent_gen("Heat transfer") + ...
    case_data.ent_gen("Species diffusion") + ...
    case_data.ent_gen("Flow viscosity")));
ax1  = subplot(2,1,1);
contourf(X, R, Y1, 15);
colormap(jet);

set(gcf, 'position', [100 100 2400 600]);
xlim([0 6]);
ylim([0 1]);
cb1 = colorbar('FontSize', 16);

ax = gca;
ax.YAxis.FontSize = 16;
ax.XAxis.FontSize = 16;

set(ylabel("$$r/R \ [-]$$", 'FontSize', 20), 'Interpreter', 'latex');
set(xlabel("$$x/R \ [-]$$", 'FontSize', 20), 'Interpreter', 'latex');
set(title("$$\sigma_\mathrm{ref} \ [\mathrm{W \, m^{-3} \, K^{-1}}]$$", 'FontSize', 23), 'Interpreter', 'latex');
axis equal;

case_data = load('D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\OneChanged\saved_results\065-05e2500-0900-005-020.mat');
Y2 = transpose((case_data.ent_gen("Chemical reactions") + ...
    case_data.ent_gen("Heat transfer") + ...
    case_data.ent_gen("Species diffusion") + ...
    case_data.ent_gen("Flow viscosity")));
ax2 = subplot(2,1,2);
contourf(X, R, Y2, 15);
colormap(jet);

set(gcf, 'position', [100 100 2400 600]);
xlim([0 6]);
ylim([0 1]);

ax = gca;
ax.YAxis.FontSize = 16;
ax.XAxis.FontSize = 16;

set(ylabel("$$r/R \ [-]$$", 'FontSize', 20), 'Interpreter', 'latex');
set(xlabel("$$x/R \ [-]$$", 'FontSize', 20), 'Interpreter', 'latex');
set(title("$$\sigma \ [\mathrm{W \, m^{-3} \, K^{-1}}]$$", 'FontSize', 23), 'Interpreter', 'latex');
axis equal;

allY = [Y1(:); Y2(:)];           % combine all data
clim = [min(allY), max(allY)];   % common color scale
set([ax1, ax2], 'CLim', clim);

cb1.Position = [0.88 0.136000000000000 0.02 0.765];


%% AAAAAAAAAAA
path = "D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\OneChanged\saved_results\";
ref_case = load(path + "reference_case.mat");
reference_porosity = 0.65;
reference_catalyst_density = 2.5e5;
reference_wall_temp = 900;
reference_inlet_v = 0.25;
reference_SC = 2;

porosities = linspace(0.5, 0.95, 20);
% catalyst_densities = linspace(10^3, 10^6, 20);
catalyst_densities = arrayfun(@(x) 10^3 * 10^(3 * x / 19), 0:19);
wall_temps = linspace(700, 1100, 20);
inlet_vs = linspace(0.05, 0.65, 20);
SCs = linspace(1.5, 4, 20);

n = 20;
temp_ent_gen_h2 = zeros([1, n]);
temp_mcr = zeros([1, n]);
temp_h2_prod = zeros([1, n]);
temp_ent_gen = zeros([1, n]);

for i = 1:n
    file_path = get_file_path(reference_porosity, reference_catalyst_density, reference_wall_temp, inlet_vs(i), reference_SC, path);
    S = load(file_path);

    temp_ent_gen_h2(i) = S.entropy_per_hydrogen;
    temp_mcr(i) = S.methane_conversion_rate;
    temp_h2_prod(i) = S.m_h2_out;
    temp_ent_gen(i) = S.ent_gen_total;
end

figure
hold on;
box on;
set(gcf,'position', [100 100 600 400]);
pbaspect([6 4 1]);
axis padded;
set(gca,'fontsize', 11.5);

xlabel("$$U_\mathrm{inlet} \ [\mathrm{m \, s^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
ylabel("$$y / y_\mathrm{ref} \ [-]$$", "Interpreter", "latex", "FontSize", 18);

plot(inlet_vs, temp_ent_gen_h2 / ref_case.entropy_per_hydrogen, 'Marker', '.', 'MarkerSize', 15, 'DisplayName', 'Specific entropy generation', 'Color', 'b');
p1 = plot(reference_inlet_v, 1, 'b*', 'MarkerSize', 9);

plot(inlet_vs, temp_mcr / ref_case.methane_conversion_rate, 'Marker', '.', 'MarkerSize', 15, 'DisplayName', 'Methane conversion rate', 'Color', 'r');
p2 = plot(reference_inlet_v, 1, 'r*', 'MarkerSize', 9);

plot(inlet_vs, temp_h2_prod / ref_case.m_h2_out, 'Marker', '.', 'MarkerSize', 15, 'DisplayName', 'Hydrogen production rate', 'Color', 'g');
p3 = plot(reference_inlet_v, 1, 'g*', 'MarkerSize', 9);

p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
p3.Annotation.LegendInformation.IconDisplayStyle = 'off';