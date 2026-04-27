% fun = @(x)[norm(x-[1,2])^2;norm(x+[2,1])^2];
% lb = [-4,-4];
% ub = -lb;
% options = optimoptions('paretosearch', 'Display', 'iter', 'MaxTime', 28800, 'UseParallel', true, ...
%     'PlotFcn',{'psplotparetof' 'psplotparetox'}, 'ParetoSetSize', 15, 'MaxFunctionEvaluations', 450);
% rng default % For reproducibility
% [x,fval,exitflag,output,residuals] = paretosearch(fun,2,[],[],[],[],lb,ub,[],options);

path = 'D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\ParetoSearch\res';
fun = @(x) funnn(x);
lb = [0.6, 800, 0.05, 2.0];
ub = [0.95, 1100, 0.65, 4.0];
rng default % For reproducibility
options = optimoptions('paretosearch', 'Display', 'iter', 'MaxTime', 36000, 'UseParallel', true, ...
    'PlotFcn',{'psplotparetof' 'psplotparetox'}, 'ParetoSetSize', 30);
[x,fval,exitflag,output,residuals] = paretosearch(fun,4,[],[],[],[],lb,ub,[],options);
save('D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\ParetoSearch\poc3.mat', ...
    "x", "fval", "exitflag", "output", "residuals");

function a = funnn(x)
    a = SteamReformerFunction_RET(x(1), 5.3448e6 * (1 - x(1)), x(2), x(3), x(4), path);
    a(1) = a(1) / 10 ^ 5;
    a(2) = a(2) * 10 ^ 5;
end

function y = funn(t, acc)
    porosity = t.por;
    catalyst_density = 5.3448e6 * (1 - porosity);
    temperature = t.temp;
    inlet_v = t.v;
    SC = t.SC;
    res = SteamReformerFunction_RET(porosity, catalyst_density, temperature, inlet_v, SC, path);
    y = sum(res .* acc);
end

% porosity = optimizableVariable("por", [0.6 0.8]);
% wall_temp = optimizableVariable("temp", [700 1100]);
% v_inlet = optimizableVariable("v", [0.05 0.65]);
% SC = optimizableVariable("SC", [1.5 4.0]);
% 
% res = bayesopt(@(x) funn(x, [2e-4 1e4 0]), [porosity wall_temp v_inlet SC], IsObjectiveDeterministic=true, UseParallel=true, MaxObjectiveEvaluations=40, MaxTime=2400);
% save('D:\AGH\MATLABPROJECT\V2\Collocated-CFD\Simulations\Results\ParetoSearch\bayes.mat', ...
%     "res");