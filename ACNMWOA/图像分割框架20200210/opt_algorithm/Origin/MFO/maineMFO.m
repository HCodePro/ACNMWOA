function [gBest, gbestvalue, FEcount, etime, MaxFEs, Convergence_curve, iter] = maineMFO(Lmax1, level, Pxy, MaxFEs, Population_size, objectiveFunction, Iteration)
    % d = level-1;
    d = level;
    NumberOfThresholdValues = 2 * d;
    thdim = NumberOfThresholdValues;
    Lmax = Lmax1;
    XVmax = Lmax * ones(1, NumberOfThresholdValues);
    XVmin = zeros(1, NumberOfThresholdValues);
    % Population_size=(thdim/2)*10;
    % MaxFEs=(thdim/2)*10000;
    % objectiveFunction='renyi2d';
    tic

    if MaxFEs == 0
        [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = IterationMFO(objectiveFunction, Population_size, thdim, XVmin, XVmax, MaxFEs, Pxy, Iteration);
    else
        [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = EvolutionMFO(objectiveFunction, Population_size, thdim, XVmin, XVmax, MaxFEs, Pxy, Iteration);
    end

    toc
    gBest = Lbest;
    gbestvalue = Fbest;
    FEcount = FE;
    etime = toc;
end
