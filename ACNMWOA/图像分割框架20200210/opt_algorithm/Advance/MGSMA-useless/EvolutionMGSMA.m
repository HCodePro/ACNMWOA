% function [BestSol,Convergence_curve]=DE( SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = EvolutionMGSMA(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    %% Initialize position
    dim = thdim / 2;
    Fbest = -inf;

    bestPositions = zeros(1, 2 * dim);
    Destination_fitness = -inf;
    AllFitness = -inf * ones(N, 1);
    weight = ones(N, 2 * dim);
    %% Initialize the set of random solutions
    Convergence_curve = [];
    it = 1;
    FEs = 0;
    X = random_initialization(N, dim, ub, lb);

    for i = 1:N
        Flag4ub = X(i, :) > ub;
        Flag4lb = X(i, :) < lb;
        X(i, :) = (X(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        %     AllFitness(i) = fobj(X(i,:));
        [AllFitness(i), FEs, X(i, :)] = sigle_evaluation(X(i, :), dim, thdim, fname, Pxy, FEs);
    end

    [SmellOrder, SmellIndex] = sort(AllFitness, 'descend'); %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);
    bestPositions = X(SmellIndex(1), :);
    Destination_fitness = bestFitness;
    ibslo = 0.8;
    %% Main loop
    while FEs < MaxFEs || it <= Iteration
        [X1] = func_MVO(AllFitness, X, bestPositions, N, FEs, Iteration, ub, lb);

        for i = 1:N
            Flag4ub = X1(i, :) > ub;
            Flag4lb = X1(i, :) < lb;
            X1(i, :) = (X1(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            [temp_fit, FEs, X1(i, :)] = sigle_evaluation(X1(i, :), dim, thdim, fname, Pxy, FEs);

            if temp_fit > AllFitness(i)
                AllFitness(i) = temp_fit;
                X(i, :) = X1(i, :);
            end

        end

        %% GK
        %         [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
        l = ceil(rand * N);

        for t = 1:N

            for i = 1:2 * dim
                D = 0;

                for r = 1:N
                    D = D + abs(X(l, i) - X(r, i));
                end

                sigma = ibslo * D / (N - 1);
                % 根据高斯核函数产生高斯随机值
                %             newpop(t).Position(i)=normrnd(s(l,i),sigma);
                pop(t, i) = X(l, i) + sigma * randn;

            end

        end

        for i = 1:N
            Flag4ub = pop(i, :) > ub;
            Flag4lb = pop(i, :) < lb;
            pop(i, :) = (pop(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            [temp_fit, FEs, pop(i, :)] = sigle_evaluation(pop(i, :), dim, thdim, fname, Pxy, FEs);

            if temp_fit > AllFitness(i)
                AllFitness(i) = temp_fit;
                X(i, :) = pop(i, :);
            end

        end

        [SmellOrder, SmellIndex] = sort(AllFitness, 'descend'); %Eq.(2.6)
        worstFitness = SmellOrder(N);
        bestFitness = SmellOrder(1);

        S = bestFitness - worstFitness + eps; % plus eps to avoid denominator zero

        %calculate the fitness weight of each slime mold
        for i = 1:N

            for j = 1:2 * dim

                if i <= (N / 2) %Eq.(2.5)
                    weight(SmellIndex(i), j) = 1 + rand() * log10((bestFitness - SmellOrder(i)) / (S) + 1);
                else
                    weight(SmellIndex(i), j) = 1 - rand() * log10((bestFitness - SmellOrder(i)) / (S) + 1);
                end

            end

        end

        %update the best fitness value and best position
        if bestFitness > Destination_fitness
            bestPositions = X(SmellIndex(1), :);
            Destination_fitness = bestFitness;
        end

        a = atanh(- (FEs / MaxFEs) + 1); %Eq.(2.4)
        b = 1 - FEs / MaxFEs;
        % Update the Position of search agents
        for i = 1:N

            if rand < 0.03 %Eq.(2.7)
                X(i, :) = (ub - lb) * rand + lb;
            else
                p = tanh(abs(AllFitness(i) - Destination_fitness)); %Eq.(2.2)
                vb = unifrnd(-a, a, 1, 2 * dim); %Eq.(2.3)
                vc = unifrnd(-b, b, 1, 2 * dim);

                for j = 1:2 * dim
                    r = rand();
                    A = randi([1, N]); % two positions randomly selected from population
                    B = randi([1, N]);

                    if r < p %Eq.(2.1)
                        X(i, j) = bestPositions(j) + vb(j) * (weight(i, j) * X(A, j) - X(B, j));
                    else
                        X(i, j) = vc(j) * X(i, j);
                    end

                end

            end

        end

        %%
        Convergence_curve(it) = Destination_fitness;

        if Fbest < Destination_fitness
            FE = FEs;
            iter = it;
        end

        Fbest = Destination_fitness;
        Lbest = bestPositions;
        it = it + 1;
    end

end
