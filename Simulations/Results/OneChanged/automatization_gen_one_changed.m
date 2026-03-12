path = "D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\OneChanged\saved_results\";

reference_porosity = 0.65;
reference_catalyst_density = 2.5e5;
reference_wall_temp = 900;
reference_inlet_v = 0.25;
reference_SC = 2;

porisities = linspace(0.5, 0.95, 20);
catalyst_densities = linspace(10^3, 10^6, 20);
wall_temps = linspace(700, 1100, 20);
inlet_vs = linspace(0.05, 0.65, 20);
SCs = linspace(1.5, 4, 20);

% Reference case
SteamReformerFunction(reference_porosity, reference_catalyst_density, ...
    reference_wall_temp, reference_inlet_v, reference_SC, ...
    path + "reference_case.mat");

fprintf("Finished the reference case\n");

% Porosity dependence
parfor i = 1:numel(porisities)
    file_path = get_file_path(porisities(i), reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, path)
    SteamReformerFunction(porisities(i), reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, file_path);
end

fprintf("Finished porosity dependence\n");

% Catalyst density dependence
parfor i = 1:numel(catalyst_densities)
    file_path = get_file_path(reference_porosity, catalyst_densities(i), reference_wall_temp, reference_inlet_v, reference_SC, path)
    SteamReformerFunction(reference_porosity, catalyst_densities(i), reference_wall_temp, reference_inlet_v, reference_SC, file_path);
end

fprintf("Finished catalyst density dependence\n");

% Wall temperature dependence
parfor i = 1:numel(wall_temps)
    file_path = get_file_path(reference_porosity, reference_catalyst_density, wall_temps(i), reference_inlet_v, reference_SC, path)
    SteamReformerFunction(reference_porosity, reference_catalyst_density, wall_temps(i), reference_inlet_v, reference_SC, file_path);
end

fprintf("Finished wall temperature dependence\n");

% Inlet velocity dependence
parfor i = 1:numel(inlet_vs)
    file_path = get_file_path(reference_porosity, reference_catalyst_density, reference_wall_temp, inlet_vs(i), reference_SC, path)
    SteamReformerFunction(reference_porosity, reference_catalyst_density, reference_wall_temp, inlet_vs(i), reference_SC, file_path);
end

fprintf("Finished inlet velocity dependence\n");

% Steam carbon ratio dependence
parfor i = 1:numel(SCs)
    file_path = get_file_path(reference_porosity, reference_catalyst_density, reference_wall_temp, reference_inlet_v, SCs(i), path)
    SteamReformerFunction(reference_porosity, reference_catalyst_density, reference_wall_temp, reference_inlet_v, SCs(i), file_path);
end

fprintf("Finished SC dependence\ns");