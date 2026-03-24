function [bestPositions, Convergence_curve] = MGSMA(N, Max_iter, lb, ub, dim, fobj)

    bestPositions = zeros(1, dim);
    Destination_fitness = inf; %change this to -inf for maximization problems
    AllFitness = inf * ones(N, 1); %record the fitness of all slime mold
    weight = ones(N, dim); %fitness weight of each slime mold
    %Initialize the set of random solutions
    X = initialization(N, dim, ub, lb);
    Convergence_curve = [];
    FEs = 0; %Number of current evaluations
    MaxFEs = Max_iter; %Maximum number of evaluations
    it = 1; %Number of iterations
    lb = ones(1, dim) .* lb; % lower boundary
    ub = ones(1, dim) .* ub; % upper boundary
    z = 0.03; % parameter

    for i = 1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub = X(i, :) > ub;
        Flag4lb = X(i, :) < lb;
        X(i, :) = (X(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        AllFitness(i) = fobj(X(i, :));
        FEs = FEs + 1;
    end

    [SmellOrder, SmellIndex] = sort(AllFitness); %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);
    bestPositions = X(SmellIndex(1), :);
    Destination_fitness = bestFitness;
    ibslo = 0.8;
    % Main loop
    while FEs < MaxFEs

        %sort the fitness
        for i = 1:N
            % Check if solutions go outside the search space and bring them back
            Flag4ub = X(i, :) > ub;
            Flag4lb = X(i, :) < lb;
            X(i, :) = (X(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            AllFitness(i) = fobj(X(i, :));
            FEs = FEs + 1;
        end

        [X1] = func_MVO(AllFitness, X, bestPositions, N, FEs, MaxFEs, ub, lb);

        for i = 1:N
            Flag4ub = X1(i, :) > ub;
            Flag4lb = X1(i, :) < lb;
            X1(i, :) = (X1(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            temp_fit = fobj(X1(i, :));
            FEs = FEs + 1;

            if temp_fit < AllFitness(i)
                AllFitness(i) = temp_fit;
                X(i, :) = X1(i, :);
            end

        end

        %% GK
        %         [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
        l = ceil(rand * N);

        for t = 1:N

            for i = 1:dim
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
            temp_fit = fobj(pop(i, :));
            FEs = FEs + 1;

            if temp_fit < AllFitness(i)
                AllFitness(i) = temp_fit;
                X(i, :) = pop(i, :);
            end

        end

        [SmellOrder, SmellIndex] = sort(AllFitness); %Eq.(2.6)
        worstFitness = SmellOrder(N);
        bestFitness = SmellOrder(1);

        S = bestFitness - worstFitness + eps; % plus eps to avoid denominator zero

        %calculate the fitness weight of each slime mold
        for i = 1:N

            for j = 1:dim

                if i <= (N / 2) %Eq.(2.5)
                    weight(SmellIndex(i), j) = 1 + rand() * log10((bestFitness - SmellOrder(i)) / (S) + 1);
                else
                    weight(SmellIndex(i), j) = 1 - rand() * log10((bestFitness - SmellOrder(i)) / (S) + 1);
                end

            end

        end

        %update the best fitness value and best position
        if bestFitness < Destination_fitness
            bestPositions = X(SmellIndex(1), :);
            Destination_fitness = bestFitness;
        end

        a = atanh(- (FEs / Max_iter) + 1); %Eq.(2.4)
        b = 1 - FEs / Max_iter;
        % Update the Position of search agents
        for i = 1:N

            if rand < z %Eq.(2.7)
                X(i, :) = (ub - lb) * rand + lb;
            else
                p = tanh(abs(AllFitness(i) - Destination_fitness)); %Eq.(2.2)
                vb = unifrnd(-a, a, 1, dim); %Eq.(2.3)
                vc = unifrnd(-b, b, 1, dim);

                for j = 1:dim
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

        Convergence_curve(it) = Destination_fitness;
        %     display(['Iteration ' num2str(it) ': Best fitness = ' num2str(Convergence_curve(it))]);

        it = it + 1;

    end

end
