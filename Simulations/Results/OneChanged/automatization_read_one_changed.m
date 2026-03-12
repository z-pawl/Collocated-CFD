path = "D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\OneChanged\saved_results\";
image_path = "C:\Users\zgolp\OneDrive\Pulpit\AGH\INŻYNIERKA\ZDJĘCIA\PARAMETRY\";

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

ref_case = load(path + "reference_case.mat");
generate_images(ref_case, image_path + "REFERENCE\", "reference");

%% Porowatość
n = 20;
temp_ent_gen_h2 = zeros([1, n]);
temp_mcr = zeros([1, n]);
temp_h2_prod = zeros([1, n]);
temp_ent_gen = zeros([1, n]);

for i = 1:n
    file_path = get_file_path(porosities(i), reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, path);
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

xlabel("$$\epsilon_0 \ [-]$$", "Interpreter", "latex", "FontSize", 18);

yyaxis left;
ylabel("$$s_\mathrm{gen,H_2} \ [\mathrm{J \, kg^{-1} \, K^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
plot(porosities, temp_ent_gen_h2, 'Marker', '.', 'MarkerSize', 15);
plot(reference_porosity, ref_case.entropy_per_hydrogen, 'b*', 'MarkerSize', 9);

yyaxis right;
ylabel("$$\mathit{\Phi}_\mathrm{CH_4} \ [-]$$", "Interpreter", "latex", "FontSize", 18);
plot(porosities, temp_mcr, 'Marker', '.', 'MarkerSize', 15);
plot(reference_porosity, ref_case.methane_conversion_rate, 'r*', 'MarkerSize', 9);

exportgraphics(gcf, image_path + 'porosity_change.pdf');

figure
hold on;
box on;
set(gcf,'position', [100 100 600 400]);
pbaspect([6 4 1]);
axis padded;
set(gca,'fontsize', 11.5);

xlabel("$$\epsilon_0 \ [-]$$", "Interpreter", "latex", "FontSize", 18);

ylabel("$$\dot{m}_\mathrm{H_2,out} \ [\mathrm{kg \, s^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
plot(porosities, temp_h2_prod, 'Marker', '.', 'MarkerSize', 15);
plot(reference_porosity, ref_case.m_h2_out, 'b*', 'MarkerSize', 9);

exportgraphics(gcf, image_path + 'porosity_change_h2_prod.pdf');

file_path = get_file_path(porosities(1), reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, path);
S = load(file_path);
generate_images(S, image_path + "POROSITY\", "porosity_min");
file_path = get_file_path(porosities(20), reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, path);
S = load(file_path);
generate_images(S, image_path + "POROSITY\", "porosity_max");

%% Gęstość katalizatora
n = 20;
temp_ent_gen_h2 = zeros([1, n]);
temp_mcr = zeros([1, n]);
temp_h2_prod = zeros([1, n]);
temp_ent_gen = zeros([1, n]);

for i = 1:n
    file_path = get_file_path(reference_porosity, catalyst_densities(i), reference_wall_temp, reference_inlet_v, reference_SC, path + "log_cat\");
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

xlabel("$$\dot{w}_\mathrm{cat} \ [\mathrm{g \, m^{-3}}]$$", "Interpreter", "latex", "FontSize", 18);
xscale log;

yyaxis left;
ylabel("$$s_\mathrm{gen,H_2} \ [\mathrm{J \, kg^{-1} \, K^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
plot(catalyst_densities, temp_ent_gen_h2, 'Marker', '.', 'MarkerSize', 15);
plot(reference_catalyst_density, ref_case.entropy_per_hydrogen, 'b*', 'MarkerSize', 9);

yyaxis right;
ylabel("$$\mathit{\Phi}_\mathrm{CH_4} \ [-]$$", "Interpreter", "latex", "FontSize", 18);
plot(catalyst_densities, temp_mcr, 'Marker', '.', 'MarkerSize', 15);
plot(reference_catalyst_density, ref_case.methane_conversion_rate, 'r*', 'MarkerSize', 9);

exportgraphics(gcf, image_path + 'catalyst_density_change.pdf');

figure
hold on;
box on;
set(gcf,'position', [100 100 600 400]);
pbaspect([6 4 1]);
axis padded;
set(gca,'fontsize', 11.5);

xlabel("$$\dot{w}_\mathrm{cat} \ [\mathrm{g \, m^{-3}}]$$", "Interpreter", "latex", "FontSize", 18);
xscale log;

ylabel("$$\dot{m}_\mathrm{H_2,out} \ [\mathrm{kg \, s^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
plot(catalyst_densities, temp_h2_prod, 'Marker', '.', 'MarkerSize', 15);
plot(reference_catalyst_density, ref_case.m_h2_out, 'b*', 'MarkerSize', 9);

exportgraphics(gcf, image_path + 'catalyst_density_change_h2_prod.pdf');

file_path = get_file_path(reference_porosity, catalyst_densities(1), reference_wall_temp, reference_inlet_v, reference_SC, path + "log_cat\");
S = load(file_path);
generate_images(S, image_path + "CATALYST DENSITY\", "catalyst_density_min");
file_path = get_file_path(reference_porosity, catalyst_densities(20), reference_wall_temp, reference_inlet_v, reference_SC, path + "log_cat\");
S = load(file_path);
generate_images(S, image_path + "CATALYST DENSITY\", "catalyst_density_max");

%% Temperatura ściany
n = 20;
temp_ent_gen_h2 = zeros([1, n]);
temp_mcr = zeros([1, n]);
temp_h2_prod = zeros([1, n]);
temp_ent_gen = zeros([1, n]);

for i = 1:n
    file_path = get_file_path(reference_porosity, reference_catalyst_density, wall_temps(i), reference_inlet_v, reference_SC, path);
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

xlabel("$$T_\mathrm{BD} \ [\mathrm{K}]$$", "Interpreter", "latex", "FontSize", 18);

yyaxis left;
ylabel("$$s_\mathrm{gen,H_2} \ [\mathrm{J \, kg^{-1} \, K^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
plot(wall_temps, temp_ent_gen_h2, 'Marker', '.', 'MarkerSize', 15);
plot(reference_wall_temp, ref_case.entropy_per_hydrogen, 'b*', 'MarkerSize', 9);

yyaxis right;
ylabel("$$\mathit{\Phi}_\mathrm{CH_4} \ [-]$$", "Interpreter", "latex", "FontSize", 18);
plot(wall_temps, temp_mcr, 'Marker', '.', 'MarkerSize', 15);
plot(reference_wall_temp, ref_case.methane_conversion_rate, 'r*', 'MarkerSize', 9);

exportgraphics(gcf, image_path + 'temperature_change.pdf');

figure
hold on;
box on;
set(gcf,'position', [100 100 600 400]);
pbaspect([6 4 1]);
axis padded;
set(gca,'fontsize', 11.5);

xlabel("$$T_\mathrm{BD} \ [\mathrm{K}]$$", "Interpreter", "latex", "FontSize", 18);

ylabel("$$\dot{m}_\mathrm{H_2,out} \ [\mathrm{kg \, s^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
plot(wall_temps, temp_h2_prod, 'Marker', '.', 'MarkerSize', 15);
plot(reference_wall_temp, ref_case.m_h2_out, 'b*', 'MarkerSize', 9);

exportgraphics(gcf, image_path + 'temperature_change_h2_prod.pdf');

file_path = get_file_path(reference_porosity, reference_catalyst_density, wall_temps(1), reference_inlet_v, reference_SC, path);
S = load(file_path);
generate_images(S, image_path + "TEMPERATURE\", "temperature_min");
file_path = get_file_path(reference_porosity, reference_catalyst_density, wall_temps(20), reference_inlet_v, reference_SC, path);
S = load(file_path);
generate_images(S, image_path + "TEMPERATURE\", "temperature_max");

%% Prędkość na wlocie
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

yyaxis left;
ylabel("$$s_\mathrm{gen,H_2} \ [\mathrm{J \, kg^{-1} \, K^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
plot(inlet_vs, temp_ent_gen_h2, 'Marker', '.', 'MarkerSize', 15);
plot(reference_inlet_v, ref_case.entropy_per_hydrogen, 'b*', 'MarkerSize', 9);

yyaxis right;
ylabel("$$\mathit{\Phi}_\mathrm{CH_4} \ [-]$$", "Interpreter", "latex", "FontSize", 18);
plot(inlet_vs, temp_mcr, 'Marker', '.', 'MarkerSize', 15);
plot(reference_inlet_v, ref_case.methane_conversion_rate, 'r*', 'MarkerSize', 9);

exportgraphics(gcf, image_path + 'velocity_change.pdf');

figure
hold on;
box on;
set(gcf,'position', [100 100 600 400]);
pbaspect([6 4 1]);
axis padded;
set(gca,'fontsize', 11.5);

xlabel("$$U_\mathrm{inlet} \ [\mathrm{m \, s^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);

ylabel("$$\dot{m}_\mathrm{H_2,out} \ [\mathrm{kg \, s^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
plot(inlet_vs, temp_h2_prod, 'Marker', '.', 'MarkerSize', 15);
plot(reference_inlet_v, ref_case.m_h2_out, 'b*', 'MarkerSize', 9);

exportgraphics(gcf, image_path + 'velocity_change_h2_prod.pdf');

file_path = get_file_path(reference_porosity, reference_catalyst_density, reference_wall_temp, inlet_vs(1), reference_SC, path);
S = load(file_path);
generate_images(S, image_path + "VELOCITY\", "velocity_min");
file_path = get_file_path(reference_porosity, reference_catalyst_density, reference_wall_temp, inlet_vs(20), reference_SC, path);
S = load(file_path);
generate_images(S, image_path + "VELOCITY\", "velocity_max");

%% SC
n = 20;
temp_ent_gen_h2 = zeros([1, n]);
temp_mcr = zeros([1, n]);
temp_h2_prod = zeros([1, n]);
temp_ent_gen = zeros([1, n]);

for i = 1:n
    file_path = get_file_path(reference_porosity, reference_catalyst_density, reference_wall_temp, reference_inlet_v, SCs(i), path);
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

xlabel("$$SC \ [-]$$", "Interpreter", "latex", "FontSize", 18);

yyaxis left;
ylabel("$$s_\mathrm{gen,H_2} \ [\mathrm{J \, kg^{-1} \, K^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
plot(SCs, temp_ent_gen_h2, 'Marker', '.', 'MarkerSize', 15);
plot(reference_SC, ref_case.entropy_per_hydrogen, 'b*', 'MarkerSize', 9);

yyaxis right;
ylabel("$$\mathit{\Phi}_\mathrm{CH_4} \ [-]$$", "Interpreter", "latex", "FontSize", 18);
plot(SCs, temp_mcr, 'Marker', '.', 'MarkerSize', 15);
plot(reference_SC, ref_case.methane_conversion_rate, 'r*', 'MarkerSize', 9);

exportgraphics(gcf, image_path + 'sc_change.pdf');

figure
hold on;
box on;
set(gcf,'position', [100 100 600 400]);
pbaspect([6 4 1]);
axis padded;
set(gca,'fontsize', 11.5);

xlabel("$$SC \ [-]$$", "Interpreter", "latex", "FontSize", 18);

ylabel("$$\dot{m}_\mathrm{H_2,out} \ [\mathrm{kg \, s^{-1}}]$$", "Interpreter", "latex", "FontSize", 18);
plot(SCs, temp_h2_prod, 'Marker', '.', 'MarkerSize', 15);
plot(reference_SC, ref_case.m_h2_out, 'b*', 'MarkerSize', 9);

exportgraphics(gcf, image_path + 'sc_change_h2_prod.pdf');

file_path = get_file_path(reference_porosity, reference_catalyst_density, reference_wall_temp, reference_inlet_v, SCs(1), path);
S = load(file_path);
generate_images(S, image_path + "SC\", "SC_min");
file_path = get_file_path(reference_porosity, reference_catalyst_density, reference_wall_temp, reference_inlet_v, SCs(20), path);
S = load(file_path);
generate_images(S, image_path + "SC\", "SC_max");