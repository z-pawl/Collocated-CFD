fun = @(x) funnn(x);
A = [0, 0, 0, -0.4, -1;
     0, 0, 0,  0.4,  1];
b = [-750; 1100];
lb = [0.6, 0.05, 750, -875, 750];
ub = [0.95, 0.65, 1100, 875, 1100];
rng default % For reproducibility
options = optimoptions('ga', 'Display', 'iter', 'MaxTime', 36000, 'UseParallel', true, ...
    'PlotFcn',{'gaplotpareto'}, 'PopulationSize', 30);
[x,fval,exitflag,output,population,scores] = gamultiobj(fun,5,A,b,[],[],lb,ub,[],options);
save('D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\TemperatureProfileOptimizationgenopt.mat', ...
    "x", "fval", "exitflag", "output", "population", "scores");

function a = funnn(x)
    xs = 0.4 / 80 * ((0:79)' + 0.5);
    wall_temps = x(4) * xs + x(5);
    y = SteamReformer_ChangingWallTemp(x(1), x(2), x(3), wall_temps);
    a = [y(1) y(2) y(5)];
end