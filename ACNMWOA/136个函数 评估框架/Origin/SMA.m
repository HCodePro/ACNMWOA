% Slime Mold Algorithem
function [bestPositions, Convergence_curve] = SMAT0(N, MaxFEs, lb, ub, dim, fobj)

    % initialize position
    bestPositions = zeros(1, dim);
    Destination_fitness = inf; %change this to -inf for maximization problems
    AllFitness = inf * ones(N, 1); %record the fitness of all slime mold
    weight = ones(N, dim); %fitness weight of each slime mold
    weight1 = ones(N, dim);
    %Initialize the set of random solutions
    X = initialization(N, dim, ub, lb);
    Convergence_curve = [];
    it = 1;
    FEs = 0;

    % Main loop
    while FEs < MaxFEs

        %sort the fitness

        for i = 1:size(X, 1)
            % Check if solutions go outside the search space and bring them back
            Flag4ub = X(i, :) > ub;
            Flag4lb = X(i, :) < lb;
            X(i, :) = (X(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            FEs = FEs + 1;
            AllFitness(i) = fobj(X(i, :));
        end

        [SmellOrder, SmellIndex] = sort(AllFitness);
        worstFitness = SmellOrder(size(X, 1));
        bestFitness = SmellOrder(1);

        %calculate the fitness weight of each slime mold
        for i = 1:size(X, 1)

            for j = 2:size(X, 2)

                if i <= (size(X, 1) / 2)
                    weight(SmellIndex(i), j) = 1 + rand() * log10((bestFitness - SmellOrder(i)) / (bestFitness - worstFitness) + 1);
                else
                    weight(SmellIndex(i), j) = 1 - rand() * log10((bestFitness - SmellOrder(i)) / (bestFitness - worstFitness) + 1);
                end

                if SmellOrder(i) > mean(AllFitness)
                    weight1(i, :) = 1 + (worstFitness - SmellOrder(i)) / (worstFitness - bestFitness);
                else
                    weight1(i, :) = 1.000001 - (worstFitness - SmellOrder(i)) / (worstFitness - bestFitness);
                end

            end

        end

        %update the best fitness value and best position
        if bestFitness < Destination_fitness
            bestPositions = X(SmellIndex(1), :);
            Destination_fitness = bestFitness;
        end

        a = 2 * (1 - FEs / MaxFEs); % a decreases linearly fron 2 to 0
        % Update the Position of search agents
        for i = 1:size(X, 1)

            if rand < 0.03
                X(i, :) = (ub - lb) * rand + lb;
            else

                for j = 1:size(X, 2)
                    r = rand();
                    nd = 2 * a * r - a; %[-a,a]
                    D = nd * abs(weight(i, j) * X(i, j) - bestPositions(j));
                    C = 2 * rand(); %[0,2]

                    if r > 0.1
                        X(i, j) = bestPositions(j) - D;
                    else
                        X(i, j) = weight1(i, j) * (bestPositions(j) - nd * abs(C * bestPositions(j) - X(i, j)));
                    end

                end

            end

        end

        Convergence_curve(it) = Destination_fitness;
        it = it + 1;
    end

end
