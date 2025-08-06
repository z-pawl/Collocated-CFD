function cost = optimization_wrapper(x)
    % Input vector is as follows: temp_inlet, temp_wall, v_inlet, SC

    addpath(genpath('D:\AGH\MATLABPROJECT\V2\Collocated-CFD'));
    
    [entropy_generation, ch4_in, ch4_out, h2_out] = steam_reformer(x(1), x(2), x(3), x(4));

    cost = 0.6 * entropy_generation / (h2_out * 1e4) + 0.4 * (ch4_in - ch4_out) / ch4_in;
end