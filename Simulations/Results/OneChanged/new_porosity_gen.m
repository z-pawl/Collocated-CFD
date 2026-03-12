path = "D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\OneChanged\saved_results\cat_por\";

reference_porosity = 0.65;
catalyst_density_normal = 5.3448 * 10^6;
reference_wall_temp = 900;
reference_inlet_v = 0.25;
reference_SC = 2;

porisities = linspace(0.5, 0.8, 20);
wall_temps = linspace(700, 1000, 20);
inlet_vs = linspace(0.05, 0.65, 20);
SCs = linspace(1.5, 4, 20);

% Reference case
SteamReformerFunction(reference_porosity, catalyst_density_normal * (1 - reference_porosity), ...
    reference_wall_temp, reference_inlet_v, reference_SC, ...
    path + "reference_case.mat");

fprintf("Finished the reference case\n");

% Porosity dependence
parfor i = 1:numel(porisities)
    file_path = get_file_path(porisities(i), catalyst_density_normal * (1 - porisities(i)), reference_wall_temp, reference_inlet_v, reference_SC, path)
    SteamReformerFunction(porisities(i), catalyst_density_normal * (1 - porisities(i)), reference_wall_temp, reference_inlet_v, reference_SC, file_path);
end

fprintf("Finished porosity dependence\n");

% Wall temperature dependence
parfor i = 1:numel(wall_temps)
    file_path = get_file_path(reference_porosity, catalyst_density_normal * (1 - reference_porosity), wall_temps(i), reference_inlet_v, reference_SC, path)
    SteamReformerFunction(reference_porosity, catalyst_density_normal * (1 - reference_porosity), wall_temps(i), reference_inlet_v, reference_SC, file_path);
end

fprintf("Finished wall temperature dependence\n");

% Inlet velocity dependence
parfor i = 1:numel(inlet_vs)
    file_path = get_file_path(reference_porosity, catalyst_density_normal * (1 - reference_porosity), reference_wall_temp, inlet_vs(i), reference_SC, path)
    SteamReformerFunction(reference_porosity, catalyst_density_normal * (1 - reference_porosity), reference_wall_temp, inlet_vs(i), reference_SC, file_path);
end

fprintf("Finished inlet velocity dependence\n");

% Steam carbon ratio dependence
parfor i = 1:numel(SCs)
    file_path = get_file_path(reference_porosity, catalyst_density_normal * (1 - reference_porosity), reference_wall_temp, reference_inlet_v, SCs(i), path)
    SteamReformerFunction(reference_porosity, catalyst_density_normal * (1 - reference_porosity), reference_wall_temp, reference_inlet_v, SCs(i), file_path);
end

fprintf("Finished SC dependence\ns");