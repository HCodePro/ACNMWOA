function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = EHGS(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)

    if Iteration ~= 0
        MaxFEs = Iteration;
    end

    dim = thdim / 2;
    Fbest = -inf;
    % disp('HGS is now tackling your problem')
    lu = [lb .* ones(1, 2 * dim); ub .* ones(1, 2 * dim)];
    % tic
    % initialize position
    bestPositions = zeros(1, 2 * dim);
    tempPosition = zeros(N, 2 * dim);

    Destination_fitness = -inf; %change this to -inf for maximization problems
    AllFitness_GWO = -inf .* ones(1, N);
    Worstest_fitness = inf;
    AllFitness = -inf * ones(1, N); %record the fitness of all positions
    VC1 = ones(1, N); %record the variation control of all positions

    weight3 = ones(N, 2 * dim); %hungry weight of each position
    weight4 = ones(N, 2 * dim); %hungry weight of each position

    %Initialize the set of random solutions
    X = random_initialization(N, dim, ub, lb);
    Convergence_curve = [];
    it = 1; %Number of iterations
    FEs = 0;

    for i = 1:N
        % Check if solutions go outside the search space and bring them back
        X(i, :) = BC(X(i, :), lb, ub);
        %     AllFitness(i) = fobj(X(i,:));
        if FEs < MaxFEs
            [AllFitness(i), FEs, X(i, :)] = sigle_evaluation(X(i, :), dim, thdim, fname, Pxy, FEs, Iteration);
        else
            break;
        end

        % FEs = FEs + 1;
    end

    % Personal best fitness and position obtained by each wolf
    pBestScore = AllFitness;
    pBest = X;
    neighbor = zeros(N, N);

    hungry = zeros(1, size(X, 1)); %record the hungry of all positions
    count = 0;

    % Main loop
    while FEs <= MaxFEs + 2

        if Iteration ~= 0
            FEs = it;
        end

        %% CrissCross
        p2 = 0.1;
        [X, AllFitness, FEs] = CC_p(X, AllFitness, dim, lb, ub, p2, thdim, fname, Pxy, FEs, Iteration, MaxFEs);
        %%
        VC2 = 0.03; %The variable of variation control

        sumHungry = 0; %record the sum of each hungry
        [AllFitnessSorted, IndexSorted] = sort(AllFitness, 'descend');
        bestFitness = AllFitnessSorted(1);
        worstFitness = AllFitnessSorted(size(X, 1));

        %update the best fitness value and best position
        if bestFitness > Destination_fitness
            bestPositions = X(IndexSorted(1), :);
            Destination_fitness = bestFitness;
            count = 0;
        end

        if worstFitness < Worstest_fitness
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
        shrink = 2 * (1 - FEs / MaxFEs); % a decreases linearly from 2 to 0

        for i = 1:size(X, 1)

            if rand < VC2

                for j = 1:size(X, 2)
                    X_GWO(i, j) = X(i, j) * (1 + randn(1));
                end

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

        %% sort the fitness
        for i = 1:size(X, 1)
            % Check if solutions go outside the search space and bring them back
            X_GWO(i, :) = boundConstraint(X_GWO(i, :), X(i, :), lu);
            %         AllFitness_GWO(i) = fobj(X_GWO(i,:));
            if FEs < MaxFEs
                [AllFitness_GWO(i), FEs, X_GWO(i, :)] = sigle_evaluation(X_GWO(i, :), dim, thdim, fname, Pxy, FEs, Iteration);
            else
                break;
            end

            % FEs = FEs + 1;
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
                random_Idx_neighbor = randi(Idx_t, 1, 2 * dim);

                for d = 1:2 * dim
                    X_DLH(t, d) = X(t, d) + rand .* (X(Idx(random_Idx_neighbor(d)), d) - X(r1(t), d)); % Equation (12)
                end

                %             X_DLH(t,:) = boundConstraint(X_DLH(t,:), X(t,:), lu);
                X_DLH(t, :) = BC(X_DLH(t, :), lb, ub);
                %         AllFitness_DLH(t) = fobj(X_DLH(t,:));
                if FEs < MaxFEs
                    [AllFitness_DLH(t), FEs, X_DLH(t, :)] = sigle_evaluation(X_DLH(t, :), dim, thdim, fname, Pxy, FEs, Iteration);
                else
                    break;
                end

            end

        end

        %% Selection
        tmp = AllFitness_GWO > AllFitness_DLH; % Equation (13)
        tmp_rep = repmat(tmp', 1, 2 * dim);

        tmpFit = tmp .* AllFitness_GWO + (1 - tmp) .* AllFitness_DLH;
        tmpPositions = tmp_rep .* X_GWO + (1 - tmp_rep) .* X_DLH;
        %% Updating
        tmp = pBestScore >= tmpFit; % Equation (13)
        tmp_rep = repmat(tmp', 1, 2 * dim);

        pBestScore = tmp .* pBestScore + (1 - tmp) .* tmpFit;
        pBest = tmp_rep .* pBest + (1 - tmp_rep) .* tmpPositions;

        AllFitness = pBestScore;
        X = pBest;
        %%
        neighbor = zeros(N, N);

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

    % toc
end

function X = BC(X, lb, ub)
    Flag4ub = X > ub;
    Flag4lb = X < lb;
    X = (X .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
end

function [X, fitness, FEs] = CC_p(X, fitness, dim, lb, ub, p2, thdim, fname, Pxy, FEs, Iteration, MaxFEs)
    %%
    % Criss
    Mhc = zeros(size(X, 1), 2 * dim);
    Bhc = randperm(size(X, 1));

    for i = 1:(size(X, 1) / 2)
        no1 = Bhc(2 * i - 1);
        no2 = Bhc(2 * i);

        for j = 1:2 * dim
            r1 = unifrnd(0, 1); %生成服从均匀分布的0-1的随机数
            r2 = unifrnd(0, 1); %生成服从均匀分布的0-1的随机数
            c1 = (rand(1) * 2) - 1; %生成服从均匀分布的-1到1的随机数
            c2 = (rand(1) * 2) - 1;
            Mhc(no1, j) = r1 * X(no1, j) + (1 - r1) * X(no2, j) + c1 * (X(no1, j) - X(no2, j));
            Mhc(no2, j) = r2 * X(no2, j) + (1 - r2) * X(no1, j) + c2 * (X(no2, j) - X(no1, j));
        end

    end

    for i = 1:size(X, 1)
        % Check boundries
        FU = Mhc(i, :) > ub; FL = Mhc(i, :) < lb; Mhc(i, :) = (Mhc(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
        % fitness of locations
        %         fitness_mhc(i)=fobj(Mhc(i,:));
        if FEs < MaxFEs
            [fitness_mhc(i), FEs, Mhc(i, :)] = sigle_evaluation(Mhc(i, :), dim, thdim, fname, Pxy, FEs, Iteration);
        else
            break;
        end

        FEs = FEs + 1;

        if fitness(i) > fitness_mhc(i)
            X(i, :) = X(i, :);
        else
            X(i, :) = Mhc(i, :);
            fitness(i) = fitness_mhc(i);
        end

    end

    %Cross
    Bvc = randperm(2 * dim);
    Mvc = X;
    %normalization
    for i = 1:size(X, 1)
        Boundary_no = size(ub, 2); % numnber of boundaries
        % If the boundaries of all variables are equal and user enter a signle
        % number for both ub and lb
        if Boundary_no == 1
            Mvc(i, :) = (Mvc(i, :) - lb) / (ub - lb);
        end

        % If each variable has a different lb and ub
        if Boundary_no > 1

            for j = 1:2 * dim
                ub_j = ub(j);
                lb_j = lb(j);
                Mvc(i, j) = (Mvc(i, j) - lb_j) / (ub_j - lb_j);
            end

        end

    end

    %     p2 = 0.6;  %p2取0.2到0.8之间
    for i = 1:(2 * dim / 2)
        p = unifrnd(0, 1); %生成服从均匀分布的0-1的随机数

        if p < p2
            no1 = Bvc(2 * i - 1);
            no2 = Bvc(2 * i);

            for j = 1:size(X, 1)
                r = unifrnd(0, 1); %生成服从均匀分布的0-1的随机数
                Mvc(j, no1) = r * Mvc(j, no1) + (1 - r) * Mvc(j, no2);
            end

        end

    end

    for i = 1:size(X, 1)
        Boundary_no = size(ub, 2); % numnber of boundaries
        % If the boundaries of all variables are equal and user enter a signle
        % number for both ub and lb
        if Boundary_no == 1
            Mvc(i, :) = Mvc(i, :) * (ub - lb) + lb;
        end

        % If each variable has a different lb and ub
        if Boundary_no > 1

            for j = 1:2 * dim
                ub_j = ub(j);
                lb_j = lb(j);
                Mvc(i, j) = (ub_j - lb_j) * Mvc(i, j) + lb_j;
            end

        end

        % Check boundries
        FU = Mvc(i, :) > ub; FL = Mvc(i, :) < lb; Mvc(i, :) = (Mvc(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
        % fitness of locations
        %         fitness_mvc(i)=fobj(Mvc(i,:));
        if FEs < MaxFEs
            [fitness_mvc(i), FEs, Mvc(i, :)] = sigle_evaluation(Mvc(i, :), dim, thdim, fname, Pxy, FEs, Iteration);
        else
            break;
        end

        FEs = FEs + 1;

        if fitness(i) > fitness_mvc(i)
            X(i, :) = X(i, :);
        else
            X(i, :) = Mvc(i, :);
            fitness(i) = fitness_mvc(i);
        end

    end

    %%
end

function [Fitness, fes, X] = sigle_evaluation(X, dim, thdim, fname, Pxy, fes, Iteration)
    %SIGLE_CALCULATE 此处显示有关此函数的摘要
    %   此处显示详细说明
    %         Vtemp1=V(1,1:dim);
    %         Vtemp1(1,:)=sort(Vtemp1(1,:));
    %         Vtemp2=V(1,(dim+1):thdim);
    %         Vtemp2(1,:)=sort(Vtemp2(1,:));
    %         V=[Vtemp1 Vtemp2] ;
    %         FitnessV=feval(fname,V(1,:),Pxy);
    %         fes = fes + 1;
    X = round(X);
    Xtemp1 = X(:, 1:dim);

    for si = 1:size(Xtemp1, 1)
        Xtemp1(si, :) = sort(Xtemp1(si, :));
    end

    Xtemp2 = X(:, (dim + 1):thdim);

    for si = 1:size(Xtemp2, 1)
        Xtemp2(si, :) = sort(Xtemp2(si, :));
    end

    X = [Xtemp1 Xtemp2];
    Fitness = feval(fname, X(1, :), Pxy);

    if Iteration == 0
        fes = fes + 1;
    end

end
