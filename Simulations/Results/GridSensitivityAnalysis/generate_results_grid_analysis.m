path = "D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\GridSensitivityAnalysis\saved_results\";

reference_porosity = 0.65;
reference_catalyst_density = 2.5e5;
reference_wall_temp = 900;
reference_inlet_v = 0.25;
reference_SC = 2;

grid_sizes = (2:12) / 2;

parfor i = 1:numel(grid_sizes)
    file_path = get_file_path_grid_analysis(grid_sizes(i), reference_porosity, reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, path)
    SteamReformerFunction(reference_porosity, reference_catalyst_density, reference_wall_temp, reference_inlet_v, reference_SC, file_path, grid_sizes(i));
end