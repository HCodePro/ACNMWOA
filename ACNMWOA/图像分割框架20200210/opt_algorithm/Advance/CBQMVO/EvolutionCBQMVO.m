function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = EvolutionCBQMVO(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    dim = thdim / 2;

    lu = [lb; ub];
    Cn = 2 * dim;
    rand('seed', sum(100 * clock));
    CN = Cn;
    maxeval = 300000; % stop criteria
    xmeanw = (lu(1, :) + rand(1, CN) .* (lu(2, :) - lu(1, :)))'; % object parameter start point
    sigma = 0.25; minsigma = 1e-15; maxsigma = max(lu(2, :) - lu(1, :)) / sqrt(Cn); % initial step size, minimal step size
    flginiphase = 1; % Initial phase
    celambda = 4 + floor(3 * log(CN));
    mu = floor(celambda / 2);
    arweights = log((celambda + 1) / 2) - log(1:mu)';
    cc = 4 / (CN + 4); ccov = 2 / (CN + 2 ^ 0.5) ^ 2;
    cs = 4 / (CN + 4); damp = (1 - min(0.7, CN * celambda / maxeval)) / cs + 1;
    B = eye(CN); D = eye(CN); BD = B * D; C = BD * transpose(BD);
    pc = zeros(CN, 1); ps = zeros(CN, 1);
    cw = sum(arweights) / norm(arweights); chiN = CN ^ 0.5 * (1 - 1 / (4 * CN) + 1 / (21 * CN ^ 2));
    counteval = 0; flgstop = 0;
    % Boundary
    Clb = (ones(celambda, 1) * lu(1, :))';
    Cub = (ones(celambda, 1) * lu(2, :))';

    %Two variables for saving the position and inflation rate (fitness) of the best universe

    Best_universe = zeros(1, 2 * dim);
    Fbest = -inf;
    Best_universe_Inflation_rate = -inf;

    %Initialize the positions of universes
    Universes = random_initialization(N, dim, ub, lb);

    %Minimum and maximum of Wormhole Existence Probability (min and max in
    % Eq.(3.3) in the paper
    WEP_Max = 1;
    WEP_Min = 0.2;
    FEs = 0;
    Convergence_curve = [];

    %Iteration(time) counter
    Time = 1;

    pBest = Universes;
    val_pBest = -inf(1, N);
    popsize = N;
    I = 1; % max immigration rate
    E = 1; % max emigration rate
    MigrateModel = 5;
    migration_models;
    GEN = 1;
    maxGEN = MaxFEs / popsize;
    iwt = 0.9 - (1:maxGEN) * (0.7 / maxGEN);

    %Main loop
    while FEs < MaxFEs || Time <= Iteration

        %Eq. (3.3) in the paper
        WEP = WEP_Min + FEs * ((WEP_Max - WEP_Min) / MaxFEs);

        %Travelling Distance Rate (Formula): Eq. (3.4) in the paper
        TDR = 1 - ((FEs) ^ (1/6) / (MaxFEs) ^ (1/6));

        %Inflation rates (I) (fitness values)
        Inflation_rates = zeros(1, size(Universes, 1));

        for i = 1:size(Universes, 1)

            %Boundary checking (to bring back the universes inside search
            % space if they go beyoud the boundaries
            Flag4ub = Universes(i, :) > ub;
            Flag4lb = Universes(i, :) < lb;
            Universes(i, :) = (Universes(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            [fitness, FEs, Universes(i, :)] = sigle_evaluation(Universes(i, :), dim, thdim, fname, Pxy, FEs);
            Inflation_rates(1, i) = fitness;
            %Elitism
            if Inflation_rates(1, i) > Best_universe_Inflation_rate
                Best_universe_Inflation_rate = Inflation_rates(1, i);
                Best_universe = Universes(i, :);
            end

            if Inflation_rates(1, i) > val_pBest(1, i)
                val_pBest(1, i) = Inflation_rates(1, i);
                pBest(i, :) = Universes(i, :);
            end

            pBest_ind(i, :) = LearnIndex_B(val_pBest, N, 2 * dim, i, mumu, lambda);
            %  Biogeography-based Learning Strategy
            for j = 1:2 * dim
                pBest_f(i, j) = pBest(pBest_ind(i, j), j);
            end

            Flag4ub = pBest_f(i, :) > ub;
            Flag4lb = pBest_f(i, :) < lb;
            pBest_f(i, :) = (pBest_f(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            [fitness, FEs, pBest_f(i, :)] = sigle_evaluation(pBest_f(i, :), dim, thdim, fname, Pxy, FEs);
            temp = fitness;
            %Elitism
            if temp > Inflation_rates(1, i)
                Inflation_rates(1, i) = temp;
                Universes(i, :) = pBest_f(i, :);
            end

            if Inflation_rates(1, i) > Best_universe_Inflation_rate
                Best_universe_Inflation_rate = Inflation_rates(1, i);
                Best_universe = Universes(i, :);
            end

        end

        % hanyan
        for i = 1:N
            temp = ub + lb - Universes(i, :);
            mid = (ub + lb) / 2;
            dis1 = pdist2(temp, Best_universe, 'Euclidean');
            dis2 = pdist2(Universes(i, :), Best_universe, 'Euclidean');

            if dis1 > dis2
                temp1 = Universes(i, :) + rand * (mid - Universes(i, :));
            else
                temp1 = temp + rand * (mid - temp);
            end

            Flag4ub = temp1 > ub;
            Flag4lb = temp1 < lb;
            temp1 = (temp1 .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            [fitness, FEs, temp1] = sigle_evaluation(temp1, dim, thdim, fname, Pxy, FEs);

            if Inflation_rates(i) < fitness
                Universes(i, :) = temp1;
            end

        end

        [sorted_Inflation_rates, sorted_indexes] = sort(Inflation_rates);

        for newindex = 1:N
            Sorted_universes(newindex, :) = Universes(sorted_indexes(newindex), :);
        end

        %Normaized inflation rates (NI in Eq. (3.1) in the paper)
        normalized_sorted_Inflation_rates = normr(sorted_Inflation_rates);

        Universes(1, :) = Sorted_universes(N, :);

        %Update the Position of universes
        for i = 1:size(Universes, 1) %Starting from 2 since the firt one is the elite
            Back_hole_index = i;

            for j = 1:size(Universes, 2)
                r1 = rand();

                if r1 < normalized_sorted_Inflation_rates(i)
                    White_hole_index = RouletteWheelSelection(-sorted_Inflation_rates); % for maximization problem -sorted_Inflation_rates should be written as sorted_Inflation_rates

                    if White_hole_index == -1
                        White_hole_index = 1;
                    end

                    %Eq. (3.1) in the paper
                    Universes(Back_hole_index, j) = Sorted_universes(White_hole_index, j);
                end

                if (size(lb, 2) == 1)
                    %Eq. (3.2) in the paper if the boundaries are all the same
                    r2 = rand();

                    if r2 < WEP
                        r3 = rand();

                        if r3 < 0.5
                            Universes(i, j) = Best_universe(1, j) + TDR * ((ub - lb) * rand + lb);
                        end

                        if r3 > 0.5
                            Universes(i, j) = Best_universe(1, j) - TDR * ((ub - lb) * rand + lb);
                        end

                    end

                end

                if (size(lb, 2) ~= 1)
                    %Eq. (3.2) in the paper if the upper and lower bounds are
                    %different for each variables
                    r2 = rand();

                    if r2 < WEP
                        r3 = rand();

                        if r3 < 0.5
                            Universes(i, j) = Best_universe(1, j) + TDR * ((ub(j) - lb(j)) * rand + lb(j));
                        end

                        if r3 > 0.5
                            Universes(i, j) = Best_universe(1, j) - TDR * ((ub(j) - lb(j)) * rand + lb(j));
                        end

                    end

                end

            end

            % Generate and evaluate celambda offspring
            arz = randn(CN, celambda);
            arx = xmeanw * ones(1, celambda) + sigma * (BD * arz);

            I = find(arx > Cub);

            arx(I) = 2 * Cub(I) - arx(I);
            aa = find(arx(I) < Clb(I));
            arx(I(aa)) = Clb(I(aa));
            I = find(arx < Clb);
            arx(I) = 2 * Clb(I) - arx(I);
            aa = find(arx(I) > Cub(I));
            arx(I(aa)) = Cub(I(aa));
            U = arx';
            arfitness = zeros(celambda, 1);

            for i = 1:celambda
                Flag4ub = U(i, :) > ub;
                Flag4lb = U(i, :) < lb;
                U(i, :) = (U(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
                [fitness, FEs, U(i, :)] = sigle_evaluation(U(i, :), dim, thdim, fname, Pxy, FEs);
                arfitness(i) = fitness;

            end

            counteval = counteval + celambda;
            % Sort by fitness and compute weighted mean in xmeanw
            [arfitness, arindex] = sort(arfitness, 'descend'); % minimization
            xold = xmeanw; % for speed up of Eq. (14)
            xmeanw = arx(:, arindex(1:mu)) * arweights / sum(arweights);
            zmeanw = arz(:, arindex(1:mu)) * arweights / sum(arweights);

            % Adapt covariance matrix
            pc = (1 - cc) * pc + (sqrt(cc * (2 - cc)) * cw / sigma) * (xmeanw - xold); % Eq. (14)

            if ~flginiphase % do not adapt in the initial phase
                C = (1 - ccov) * C + ccov * pc * transpose(pc); % Eq. (15)
            end

            % adapt sigma
            ps = (1 - cs) * ps + (sqrt(cs * (2 - cs)) * cw) * (B * zmeanw); % Eq. (16)
            sigma = sigma * exp((norm(ps) - chiN) / chiN / damp); % Eq. (17)

            % Update B and D from C
            if mod(counteval / celambda, 1 / ccov / CN / 5) < 1
                C = triu(C) + transpose(triu(C, 1)); % enforce symmetry
                [B, D] = eig(C);
                % limit condition of C to 1e14 + 1
                if max(diag(D)) > 1e14 * min(diag(D))
                    tmp = max(diag(D)) / 1e14 - min(diag(D));
                    C = C + tmp * eye(CN); D = D + tmp * eye(CN);
                end

                D = diag(sqrt(diag(D))); % D contains standard deviations now
                BD = B * D; % for speed up only
            end % if mod

            % Adjust minimal step size
            if sigma * min(diag(D)) < minsigma ...
                    | arfitness(1) == arfitness(min(mu + 1, celambda)) ...
                    | xmeanw == xmeanw ...
                    + 0.2 * sigma * BD(:, 1 + floor(mod(counteval / celambda, CN)))
                sigma = 1.4 * sigma;

                % flgstop = 1;
            end

            if sigma > maxsigma
                sigma = maxsigma;
            end

            % Test for end of initial phase
            if flginiphase && counteval / celambda > 2 / cs

                if (norm(ps) - chiN) / chiN < 0.05 % step size is not much too small
                    flginiphase = 0;
                end

            end

            if (arfitness(1) > Best_universe_Inflation_rate)
                Best_universe_Inflation_rate = arfitness(1);
                Best_universe = U(arindex(1), :);
            end

        end

        Convergence_curve(Time) = Best_universe_Inflation_rate;

        if Fbest < Best_universe_Inflation_rate
            FE = FEs;
            iter = Time;
        end

        Fbest = Best_universe_Inflation_rate;
        Lbest = Best_universe;
        %Print the best universe details after every 50 iterations
        %     if mod(Time,50)==0
        %         display(['At iteration ', num2str(Time), ' the best universes fitness is ', num2str(Best_universe_Inflation_rate)]);
        %     end
        Time = Time + 1;
    end

end

function pBest_ind = LearnIndex_B(pfit, popsize, D, i, mu, lambda)
    % Biogeography-based exemplar generation method

    [~, sort_index] = sort(pfit, 'descend');

    for k = 1:length(sort_index)
        index(sort_index(k)) = k;
    end

    mu1 = mu(index);
    lambda1 = lambda(index);

    for k = 1:D

        if rand < lambda1(i) %  Should we immigrate?
            % Yes - Pick a solution from which to emigrate (roulette wheel selection)
            RandomNum = rand * sum(mu1);
            Select = mu1(1);
            SelectIndex = 1;

            while (RandomNum > Select) && (SelectIndex < popsize)
                SelectIndex = SelectIndex + 1;
                Select = Select + mu1(SelectIndex);
            end

            pBest_ind(k) = SelectIndex;
        else
            pBest_ind(k) = i;
        end

    end

    if all(pBest_ind == i)
        ind = randi(popsize);
        while ind == i, ind = randi(popsize); end
        pBest_ind(randi(D)) = ind;
    end

end
