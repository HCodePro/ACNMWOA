% function [Destination_fitness,bestPositions,Convergence_curve]=HGS(N,Max_iter,lb,ub,dim,fobj)
function [bestPositions, Convergence_curve] = HGS_DLH(N, MaxFEs, lb, ub, dim, fobj)
    % disp('HGS is now tackling your problem')
    lu = [lb .* ones(1, dim); ub .* ones(1, dim)];
    p2 = 0.1;
    % tic
    % initialize position
    bestPositions = zeros(1, dim);
    tempPosition = zeros(N, dim);

    Destination_fitness = inf; %change this to -inf for maximization problems
    Worstest_fitness = -inf;
    AllFitness = inf * ones(1, N); %record the fitness of all positions
    VC1 = ones(1, N); %record the variation control of all positions

    weight3 = ones(N, dim); %hungry weight of each position
    weight4 = ones(N, dim); %hungry weight of each position

    %Initialize the set of random solutions
    X = initialization(N, dim, ub, lb);
    Convergence_curve = [];
    it = 1; %Number of iterations
    FEs = 0;

    for i = 1:N
        % Check if solutions go outside the search space and bring them back
        X(i, :) = BC(X(i, :), lb, ub);
        AllFitness(i) = fobj(X(i, :));
        FEs = FEs + 1;
    end

    % Personal best fitness and position obtained by each wolf
    pBestScore = AllFitness;
    pBest = X;
    neighbor = zeros(N, N);

    hungry = zeros(1, size(X, 1)); %record the hungry of all positions
    count = 0;

    % Main loop
    while FEs <= MaxFEs
        %% CrissCross

        [X, AllFitness, FEs] = CC_p(X, AllFitness, dim, lb, ub, fobj, FEs, p2);
        %%
        VC2 = 0.03; %The variable of variation control

        sumHungry = 0; %record the sum of each hungry
        [AllFitnessSorted, IndexSorted] = sort(AllFitness);
        bestFitness = AllFitnessSorted(1);
        worstFitness = AllFitnessSorted(size(X, 1));

        %update the best fitness value and best position
        if bestFitness < Destination_fitness
            bestPositions = X(IndexSorted(1), :);
            Destination_fitness = bestFitness;
            count = 0;
        end

        if worstFitness > Worstest_fitness
            Worstest_fitness = worstFitness;
        end

        for i = 1:size(X, 1)
            %calculate the variation control of all positions
            VC1(i) = sech(abs(AllFitness(i) - Destination_fitness));
            %calculate the hungry of each position
            if Destination_fitness == AllFitness(i)
                hungry(1, i) = 0;
                count = count + 1;
                tempPosition(count, :) = X(i, :);
            else
                temprand = rand();
                c = (AllFitness(i) - Destination_fitness) / (Worstest_fitness - Destination_fitness) * temprand * 2 * (ub - lb);

                if c < 100
                    b = 100 * (1 + temprand);
                else
                    b = c;
                end

                hungry(1, i) = hungry(1, i) + max(b);
                sumHungry = sumHungry + hungry(1, i);
            end

        end

        %calculate the hungry weight of each position
        for i = 1:size(X, 1)

            for j = 2:size(X, 2)
                weight3(i, j) = (1 - exp(-abs(hungry(1, i) - sumHungry))) * rand() * 2;

                if rand() < VC2
                    weight4(i, j) = hungry(1, i) * size(X, 1) / sumHungry * rand();
                else
                    weight4(i, j) = 1;
                end

            end

        end

        % Update the Position of search agents
        shrink = 2 * (1 - FEs / MaxFEs); % a decreases linearly fron 2 to 0

        for i = 1:size(X, 1)

            if rand < VC2
                X_GWO(i, :) = X(i, j) * (1 + randn(1));
            else

                if count >= 1
                    A = randi([1, count]);

                    for j = 1:size(X, 2)
                        r = rand();
                        vb = 2 * shrink * r - shrink; %[-a,a]
                        % Moving based on the bestPosition
                        % The transformation range is controlled by weight3,bestPositions and X
                        if r > VC1(i)
                            X_GWO(i, j) = weight4(i, j) * tempPosition(A, j) + vb * weight3(i, j) * abs(tempPosition(A, j) - X(i, j));
                        else
                            X_GWO(i, j) = weight4(i, j) * tempPosition(A, j) - vb * weight3(i, j) * abs(tempPosition(A, j) - X(i, j));
                        end

                    end

                end

            end

        end

        %sort the fitness
        for i = 1:N
            % Check if solutions go outside the search space and bring them back
            X_GWO(i, :) = boundConstraint(X_GWO(i, :), X(i, :), lu);
            AllFitness_GWO(i) = fobj(X_GWO(i, :));
            FEs = FEs + 1;
        end

        %% Calculate the candiadate position Xi-DLH
        radius = pdist2(X, X_GWO, 'euclidean'); % Equation (10)
        dist_Position = squareform(pdist(X));
        r1 = randperm(N, N);

        for t = 1:N
            neighbor(t, :) = (dist_Position(t, :) <= radius(t, t));
            [~, Idx] = find(neighbor(t, :) == 1); % Equation (11)

            if isempty(Idx) == 1
                X_DLH(t, :) = X(t, :);
                AllFitness_DLH(t) = AllFitness(t);
            else
                Idx_t = size(Idx, 2);
                %         if  isempty(Idx_t)==1
                %             Idx_t=2;
                %         end
                random_Idx_neighbor = randi(Idx_t, 1, dim);

                for d = 1:dim
                    X_DLH(t, d) = X(t, d) + rand .* (X(Idx(random_Idx_neighbor(d)), d) ...
                        - X(r1(t), d)); % Equation (12)
                end

                X_DLH(t, :) = boundConstraint(X_DLH(t, :), X(t, :), lu);
                AllFitness_DLH(t) = fobj(X_DLH(t, :));
            end

        end

        %% Selection
        tmp = AllFitness_GWO < AllFitness_DLH; % Equation (13)
        tmp_rep = repmat(tmp', 1, dim);

        tmpFit = tmp .* AllFitness_GWO + (1 - tmp) .* AllFitness_DLH;
        tmpPositions = tmp_rep .* X_GWO + (1 - tmp_rep) .* X_DLH;
        %% Updating
        tmp = pBestScore <= tmpFit; % Equation (13)
        tmp_rep = repmat(tmp', 1, dim);

        pBestScore = tmp .* pBestScore + (1 - tmp) .* tmpFit;
        pBest = tmp_rep .* pBest + (1 - tmp_rep) .* tmpPositions;

        AllFitness = pBestScore;
        X = pBest;
        %%
        neighbor = zeros(N, N);

        %%
        Convergence_curve(it) = Destination_fitness;
        it = it + 1;
    end

    % toc
end

function X = BC(X, lb, ub)
    Flag4ub = X > ub;
    Flag4lb = X < lb;
    X = (X .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
end
