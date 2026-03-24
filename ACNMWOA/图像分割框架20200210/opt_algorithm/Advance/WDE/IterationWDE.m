function [Fbest, Lbest, FE, MaxEpk, Convergence_curve, iter] = IterationWDE(fname, N, thdim, low, up, MaxEpk, Pxy, Iteration)
    D = thdim / 2;
    Fbest = -inf;
    Convergence_curve = [];
    FEs = 0;
    %INITIALIZATION
    if numel(low) == 1, low = low * ones(1, 2 * D); up = up * ones(1, 2 * D); end % this line must be adapted to your problem
    % P = GenP(2*N,2*D,low,up); % see Eq.1 in [1]
    P = random_initialization(2 * N, D, low, up);
    it = 1;

    for i = 1:2 * N
        % 	fitP(i) = fobj(P(i,:));
        %     FEs=FEs+1;
        if it < Iteration
            [fitP(i), FEs, P(i, :)] = sigle_evaluation(P(i, :), D, thdim, fname, Pxy, FEs);
        else
            break
        end

    end

    % ------------------------------------------------------------------------------------------

    while FEs < MaxEpk || it <= Iteration

        j = randperm(2 * N);
        k = j(1:N);
        l = j(N + 1:2 * N);
        trialP = P(k, :);
        fitTrialP = fitP(k);

        temp = trialP; % memory

        for index = 1:N
            w = rand(N, 1) .^ 3;
            w = w ./ sum(w);
            temp(index, :) = sum(w .* P(l, :));
        end

        while 1, m = randperm(N); if sum(1:N == m) == 0, break; end, end %  bijective-morphism

        E = temp - trialP(m, :);

        %  recombination
        M = GenM(N, 2 * D);

        if rand < rand % F is dominant on problem solving success of WDE.
            % F = randn.^3;
            F = randn(1, 2 * D) .^ 3; % default
            % F = 3*randn;
            % F = (rand-.5)/rand;
            % F = 1 / gamrnd(1,0.5); % pseudo-stable walk (levy-like, simulates inverse gamma distribution; levy-distiribution)
            % F = 1 / normrnd(0,5);  % pseudo-stable walk (levy-like)
        else
            F = randn(N, 1) .^ 3;
            % F = 3*randn(N,1);
            % F = (rand(N,1)-.5)/rand;
            % F = 1 / gamrnd(1,0.5); % pseudo-stable walk (levy-like, simulates inverse gamma distribution; levy-distiribution)
            % F = 1 / normrnd(0,5);  % pseudo-stable walk (levy-like)
        end

        Trial = trialP + F .* M .* E; % re-scaling and shift

        Trial = BoundaryControl(Trial, low, up);

        for i = 1:N
            % 		fitT(i) = fobj(Trial(i,:));
            %         FEs=FEs+1;
            [fitT(i), FEs, Trial(i, :)] = sigle_evaluation(Trial(i, :), D, thdim, fname, Pxy, FEs);
        end

        ind = fitT > fitTrialP;

        trialP(ind, :) = Trial(ind, :);
        fitTrialP(ind) = fitT(ind);

        fitP(k) = fitTrialP;
        P(k, :) = trialP;

        % keep the solutions
        [bestsol, ind] = max(fitP);
        best = P(ind, :);
        Convergence_curve(it) = bestsol;

        if Fbest < bestsol
            FE = FEs;
            iter = it;
        end

        Fbest = bestsol;
        Lbest = best;
        it = it + 1;

    end %epk

    return

    function M = GenM(N, D)
        M = zeros(N, D);

        for i = 1:N
            if rand < rand, k = rand ^ 3; else, k = 1 - rand ^ 3; end
            V = randperm(D);
            j = V(1:ceil(k * D));
            M(i, j) = 1;
        end

        function pop = GenP(N, D, low, up)
            pop = ones(N, D);

            for i = 1:N

                for j = 1:D
                    pop(i, j) = rand * (up(j) - low(j)) + low(j);
                end

            end

            return

            function pop = BoundaryControl(pop, low, up)
                [popsize, dim] = size(pop);

                for i = 1:popsize

                    for j = 1:dim
                        F = rand .^ 3;
                        if pop(i, j) < low(j), pop(i, j) = low(j) + F .* (up(j) - low(j)); end
                        if pop(i, j) > up(j), pop(i, j) = up(j) + F .* (low(j) - up(j)); end
                    end

                end

                return
