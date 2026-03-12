function generate_images(case_data, image_path, name_prefix)
    radius = 0.05;
    length = 0.3;
    [X, R] = meshgrid(case_data.grid.cent_pos_x, case_data.grid.cent_pos_r);
    X = (X - X(1)) / (X(end) - X(1)) * length / radius;
    R = (R - R(1)) / (R(end) - R(1));

    %% Total generation
    figure;
    contourf(X, R, transpose((case_data.ent_gen("Chemical reactions") + ...
        case_data.ent_gen("Heat transfer") + ...
        case_data.ent_gen("Species diffusion") + ...
        case_data.ent_gen("Flow viscosity"))), 15);
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
    set(title("$$\sigma \ [\mathrm{W \, m^{-3} \, K^{-1}}]$$", 'FontSize', 23), 'Interpreter', 'latex');
    axis equal;
    
    exportgraphics(gcf, image_path + name_prefix + '_total_contour.pdf');

    %% Chemical reactions
    figure;
    contourf(X, R, transpose((case_data.ent_gen("Chemical reactions"))), 8);
    colormap(jet);
    
    set(gcf, 'position', [100 100 1800 420]);
    xlim([0 6]);
    ylim([0 1]);
    cb = colorbar('FontSize', 25);
    
    ax = gca;
    ax.YAxis.FontSize = 32;
    ax.XAxis.FontSize = 32;
    
    set(ylabel("$$r/R \ [-]$$", 'FontSize', 38), 'Interpreter', 'latex');
    set(xlabel("$$x/R \ [-]$$", 'FontSize', 38), 'Interpreter', 'latex');
    set(title("$$\sigma_\mathrm{CR} \ [\mathrm{W \, m^{-3} \, K^{-1}}]$$", 'FontSize', 46), 'Interpreter', 'latex');
    axis equal;
    
    exportgraphics(gcf, image_path + name_prefix + '_CR_contour.pdf');

    %% Heat transfer
    figure;
    contourf(X, R, transpose((case_data.ent_gen("Heat transfer"))), 8);
    colormap(jet);
    
    set(gcf, 'position', [100 100 1800 420]);
    xlim([0 6]);
    ylim([0 1]);
    cb = colorbar('FontSize', 25);
    
    ax = gca;
    ax.YAxis.FontSize = 32;
    ax.XAxis.FontSize = 32;
    
    set(ylabel("$$r/R \ [-]$$", 'FontSize', 38), 'Interpreter', 'latex');
    set(xlabel("$$x/R \ [-]$$", 'FontSize', 38), 'Interpreter', 'latex');
    set(title("$$\sigma_\mathrm{HT} \ [\mathrm{W \, m^{-3} \, K^{-1}}]$$", 'FontSize', 46), 'Interpreter', 'latex');
    axis equal;
    
    exportgraphics(gcf, image_path + name_prefix + '_HT_contour.pdf');

    %% Flow viscosity
    figure;
    contourf(X, R, transpose((case_data.ent_gen("Flow viscosity"))), 8);
    colormap(jet);
    
    set(gcf, 'position', [100 100 1800 420]);
    xlim([0 6]);
    ylim([0 1]);
    cb = colorbar('FontSize', 25);
    
    ax = gca;
    ax.YAxis.FontSize = 32;
    ax.XAxis.FontSize = 32;
    
    set(ylabel("$$r/R \ [-]$$", 'FontSize', 38), 'Interpreter', 'latex');
    set(xlabel("$$x/R \ [-]$$", 'FontSize', 38), 'Interpreter', 'latex');
    set(title("$$\sigma_\mathrm{FV} \ [\mathrm{W \, m^{-3} \, K^{-1}}]$$", 'FontSize', 46), 'Interpreter', 'latex');
    axis equal;
    
    exportgraphics(gcf, image_path + name_prefix + '_FV_contour.pdf');

    %% Species diffusion
    figure;
    contourf(X, R, transpose((case_data.ent_gen("Species diffusion"))), 8);
    colormap(jet);
    
    set(gcf, 'position', [100 100 1800 420]);
    xlim([0 6]);
    ylim([0 1]);
    cb = colorbar('FontSize', 25);
    
    ax = gca;
    ax.YAxis.FontSize = 32;
    ax.XAxis.FontSize = 32;
    
    set(ylabel("$$r/R \ [-]$$", 'FontSize', 38), 'Interpreter', 'latex');
    set(xlabel("$$x/R \ [-]$$", 'FontSize', 38), 'Interpreter', 'latex');
    set(title("$$\sigma_\mathrm{SD} \ [\mathrm{W \, m^{-3} \, K^{-1}}]$$", 'FontSize', 46), 'Interpreter', 'latex');
    axis equal;
    
    exportgraphics(gcf, image_path + name_prefix + '_SD_contour.pdf');
end