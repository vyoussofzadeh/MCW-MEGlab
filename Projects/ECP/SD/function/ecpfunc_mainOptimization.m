function ecpfunc_mainOptimization()
    bestConcordance = -inf;
    bestMethod = '';
    bestBounds = [];

    function [cost, method] = objectiveFunction(x)
        cfg.lowerBound = x(1);
        cfg.upperBound = x(2);
        [concordance, method] = ecpfunc_optimalwindows_dics(cfg);
        cost = -concordance;  % Since fmincon minimizes
        if -cost > bestConcordance
            bestConcordance = -cost;
            bestBounds = x;
            bestMethod = method;
        end
    end

    % Run the optimizer
    [bounds, fval] = fmincon(@objectiveFunction, x0, [], [], [], [], lb, ub, [], options);
    % Use the stored best values directly, avoiding re-evaluation
    fprintf('Best Method: %s, Best Bounds: [%f, %f], Max Concordance: %f\n', bestMethod, bestBounds(1), bestBounds(2), bestConcordance);
end