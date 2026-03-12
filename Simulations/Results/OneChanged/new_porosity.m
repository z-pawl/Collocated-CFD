path = "D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\OneChanged\saved_results\cat_por\";
image_path = "C:\Users\zgolp\OneDrive\Pulpit\AGH\INŻYNIERKA\ZDJĘCIA\PARAMETRY\";

reference_porosity = 0.65;
catalyst_density_normal = 5.3448 * 10^6;
reference_wall_temp = 900;
reference_inlet_v = 0.25;
reference_SC = 2;

porosities = linspace(0.5, 0.8, 20);

% Porowatość
n = 20;
temp_ent_gen_h2 = zeros([1, n]);
temp_mcr = zeros([1, n]);

for i = 1:n
    file_path = get_file_path(porosities(i), catalyst_density_normal * (1 - porosities(i)), reference_wall_temp, reference_inlet_v, reference_SC, path);
    S = load(file_path);
    fix_res(S, file_path)


    temp_ent_gen_h2(i) = S.entropy_per_hydrogen;
    temp_mcr(i) = S.methane_conversion_rate;
end

figure
hold on;
box on;
set(gcf,'position', [100 100 600 400]);
pbaspect([6 4 1]);
axis padded;

xlabel("$$\epsilon_0 \ [-]$$", "Interpreter", "latex");

yyaxis left;
ylabel("$$s_{gen,H_2} \ [J kg^{-1} K^{-1}]$$", "Interpreter", "latex");
plot(porosities, temp_ent_gen_h2, 'Marker', '.', 'MarkerSize', 15);
plot(reference_porosity, ref_case.entropy_per_hydrogen, 'b*', 'MarkerSize', 9);

yyaxis right;
ylabel("$$\Phi_{CH_4} \ [-]$$", "Interpreter", "latex");
plot(porosities, temp_mcr, 'Marker', '.', 'MarkerSize', 15);
plot(reference_porosity, ref_case.methane_conversion_rate, 'r*', 'MarkerSize', 9);

% exportgraphics(gcf, image_path + 'porosity_change.pdf');
% 
% file_path = get_file_path(porosities(1), reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, path);
% S = load(file_path);
% generate_images(S, image_path + "POROSITY\", "porosity_min");
% file_path = get_file_path(porosities(20), reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, path);
% S = load(file_path);
% generate_images(S, image_path + "POROSITY\", "porosity_max");