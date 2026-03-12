n = 200;

path = "D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\AllChanged\saved_results\";
res_structs = idk(n, path);

figure(1);
hold on;
plot(1:n, arrayfun(@(x) res_structs(x).entropy_per_hydrogen, 1:n));
plot(1:n, arrayfun(@(x) res_structs(x).methane_conversion_rate, 1:n));


function res_struct = idk(n, path)
    porosities = [0.60, 0.65, 0.70];
    catalyst_densities = [2.5e3, 2.5e4, 2.5e5];
    temps_wall = [900, 950, 1000];
    vs_inlet = [0.15, 0.20, 0.25, 0.30];
    SCs = [1.5, 2.0, 2.5, 3.0];
    
    [all_porosities, all_catalyst_densities, all_temps_wall, all_vs_inlet, all_SCs] = ndgrid(porosities, catalyst_densities, temps_wall, vs_inlet, SCs);

    opt_param = linspace(0, 1, n);
    opt_val = -inf * ones(1, n);
    opt_ind = zeros(1, n);

    for i = 1:numel(all_porosities)
        name = get_file_path(all_porosities(i), all_catalyst_densities(i), all_temps_wall(i), all_vs_inlet(i), all_SCs(i), path);
        temp = load(name);

        mcr = temp.methane_conversion_rate;
        ent_per_h2 = temp.entropy_per_hydrogen;

        for j = 1:n
            test_val = opt_param(j) * mcr - (1 - opt_param(j)) * ent_per_h2;
            if test_val > opt_val(j)
                opt_val(j) = test_val;
                opt_ind(j) = i;
            end
        end
    end

    for j = 1:n
        porosity_j = all_porosities(opt_ind(j));
        catalyst_density_j = all_catalyst_densities(opt_ind(j));
        temp_wall_j = all_temps_wall(opt_ind(j));
        v_inlet_j = all_vs_inlet(opt_ind(j));
        SC_j = all_SCs(opt_ind(j));

        best_temp = load(get_file_path(porosity_j, catalyst_density_j, temp_wall_j, v_inlet_j, SC_j, path));
        best_temp.opt_param = opt_param(j);
        best_temp.porosity = porosity_j;
        best_temp.catalyst_density = catalyst_density_j;
        best_temp.temp_wall = temp_wall_j;
        best_temp.v_inlet = v_inlet_j;
        best_temp.SC = SC_j;

        res_struct(j) = best_temp;
    end
end