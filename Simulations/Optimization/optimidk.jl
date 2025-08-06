using MATLAB, Optim

# Initialize MATLAB session
mat"disp('MATLAB started from Julia')"
mat"addpath('D:/AGH/MATLABPROJECT/V2/Collocated-CFD/Simulations/Optimization')"

# Define a Julia wrapper for the MATLAB function
function matlab_wrapper(x)
    # Send x to MATLAB
    matlab_result = mxcall(:optimization_wrapper, 1, x)
    return matlab_result[1]
end

function foo(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

function optimize_idk()
    lower = [900.0, 900.0, 0.05, 1.0]
    upper = [1200.0, 1200.0, 0.4, 3.5]
    initial_x = [1000.0, 1000.0, 0.25, 2.0]
    # lower = [-10.0, -10.0]
    # upper = [10.0, 10.0]
    # initial_x = [0.0, 0.0]
    return optimize(
        matlab_wrapper,
        lower,
        upper,
        initial_x,
        Fminbox(NelderMead()),
        Optim.Options(f_calls_limit = 100, show_trace = true))
end

res = optimize_idk()
println(res)