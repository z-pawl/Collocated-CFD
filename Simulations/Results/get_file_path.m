function name = get_file_path(porosity, catalyst_density, temp_wall, v_inlet, SC, path)
    % Nazwa
    
    % Porowatość - 3 ostatnie cyfry części całkowitej
    % 100-krotności porowatości
    s_porosity = sprintf("%03d", floor(porosity * 100));
    
    % Gęstość katalizatora - 
    % Zapisuje 2 ostatnie cyfry wykładnika i 4 ostatnie cyfry części całkowitej
    % 1000-krotności mantysy dla reprezentacji naukowej
    parts = strsplit(sprintf("%e", catalyst_density), "e+");
    exponent = str2num(parts{2});
    mantissa = str2num(parts{1});
    s_catalyst_density = sprintf("%02d", exponent) + "e" + sprintf("%04d", floor(mantissa * 1000));
    
    % Temperatura ściany - 4 ostatnie cyfry częśći całkowitej
    s_temp_wall = sprintf("%04d", floor(temp_wall));
    
    % Prędkość na wlocie - 3 ostatnie cyfry części całkowitej
    % 100-krotności prędkości na wlocie
    s_v_inlet = sprintf("%03d", floor(v_inlet * 100));
    
    % SC - 3 ostatnie cyfry części całkowitej 10-krotności SC
    s_SC = sprintf("%03d", floor(SC * 10));
    
    name = sprintf("%s-%s-%s-%s-%s.mat", [s_porosity, s_catalyst_density, s_temp_wall, s_v_inlet, s_SC]);
    name = path + name;
end