function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = IterationSMA(fname, SearchAgents_no, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    % clear
    % load('figure.mat')

    dim = thdim / 2;
    Fbest = -inf;
    N = SearchAgents_no;

    % disp('SMA_FEs is now tackling your problem')

    % initialize position
    bestPositions = zeros(1, 2 * dim);
    Destination_fitness = -inf; %change this to -inf for maximization problems
    AllFitness = -inf * ones(N, 1); %record the fitness of all slime mold
    weight = ones(N, 2 * dim); %fitness weight of each slime mold
    %Initialize the set of random solutions
    % X=initialization(N,2*dim,ub,lb);
    X = random_initialization(N, dim, ub, lb);

    Convergence_curve = [];
    FEs = 0; %Number of current evaluations
    it = 1; %Number of iterations
    lb = ones(1, 2 * dim) .* lb; % lower boundary
    ub = ones(1, 2 * dim) .* ub; % upper boundary
    z = 0.03; % parameter

    % Main loop
    while FEs < MaxFEs || it <= Iteration
        %sort the fitness
        for i = 1:N
            % Check if solutions go outside the search space and bring them back
            Flag4ub = X(i, :) > ub;
            Flag4lb = X(i, :) < lb;
            X(i, :) = (X(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            %         AllFitness(i) = fobj(X(i,:));
            %         FEs = FEs+1;
            [AllFitness(i), FEs, X(i, :)] = sigle_evaluation(X(i, :), dim, thdim, fname, Pxy, FEs, Iteration);

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

        a = atanh(- (it / Iteration) + 1); %Eq.(2.4)
        b = 1 - it / Iteration;
        % Update the Position of search agents
        for i = 1:N

            if rand < z %Eq.(2.7)
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
