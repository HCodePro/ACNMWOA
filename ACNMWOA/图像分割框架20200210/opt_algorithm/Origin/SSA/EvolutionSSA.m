% function [FoodFitness,FoodPosition,Convergence_curve]=SSA(N,Max_iter,lb,ub,dim,fobj)
function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = EvolutionSSA(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    % FoodFitness=-inf;
    %%反向学习+DE交叉变异
    %%引自Y. Wang, Z. Cai, Q. Zhang, Differential Evolution With Composite Trial Vector Generation Strategies and Control Parameters,
    %%IEEE Transactions on Evolutionary Computation, 15 (2011) 55-66.
    Fbest = -inf;
    dim = thdim / 2;
    FEs = 0;
    %Initialize the positions of salps
    SalpPositions = random_initialization(N, dim, ub, lb);
    FoodPosition = zeros(1, 2 * dim);
    FoodFitness = -inf;
    %calculate the fitness of initial salps

    for i = 1:size(SalpPositions, 1)

        if FEs < MaxFEs
            [SalpFitness(1, i), FEs, SalpPositions(i, :)] = sigle_evaluation(SalpPositions(i, :), dim, thdim, fname, Pxy, FEs);
        else
            break;
        end

    end

    [sorted_salps_fitness, sorted_indexes] = sort(SalpFitness, 'descend');

    for newindex = 1:size(SalpPositions, 1)
        Sorted_salps(newindex, :) = SalpPositions(sorted_indexes(newindex), :);
    end

    FoodPosition = Sorted_salps(1, :);
    FoodFitness = sorted_salps_fitness(1);
    Convergence_curve(1) = FoodFitness;
    it = 2;
    %Main loop
    while FEs < MaxFEs || it <= Iteration
        c1 = 2 * exp(- (4 * FEs / MaxFEs) ^ 2); % Eq. (3.2) in the paper

        for i = 1:size(SalpPositions, 1)

            if i <= N / 2

                for j = 1:size(SalpPositions, 2)
                    c2 = rand();
                    c3 = rand();
                    %%%%%%%%%%%%% % Eq. (3.1) in the paper %%%%%%%%%%%%%%
                    if c3 < 0.5
                        SalpPositions(i, j) = FoodPosition(j) + c1 * ((ub(j) - lb(j)) * c2 + lb(j));
                    else
                        SalpPositions(i, j) = FoodPosition(j) - c1 * ((ub(j) - lb(j)) * c2 + lb(j));
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end

            else
                point1 = SalpPositions(i - 1, :);
                point2 = SalpPositions(i, :);

                SalpPositions(i, :) = (point2 + point1) / 2; % % Eq. (3.4) in the paper
            end

        end

        for i = 1:size(SalpPositions, 1)
            Tp = SalpPositions(i, :) > ub;
            Tm = SalpPositions(i, :) < lb;
            SalpPositions(i, :) = (SalpPositions(i, :) .* (~(Tp + Tm))) + ub .* Tp + lb .* Tm;

            if FEs < MaxFEs
                [SalpFitness(1, i), FEs, SalpPositions(i, :)] = sigle_evaluation(SalpPositions(i, :), dim, thdim, fname, Pxy, FEs);

                if SalpFitness(1, i) > FoodFitness
                    FoodPosition = SalpPositions(i, :);
                    FoodFitness = SalpFitness(1, i);
                end

            else
                break;
            end

        end

        Convergence_curve(it) = FoodFitness;

        if Fbest < FoodFitness
            FE = FEs;
            iter = it;
        end

        Fbest = FoodFitness;
        Lbest = FoodPosition;
        it = it + 1;
        % 	FE=FEs;
    end

end

function o = Levy1(d)
    beta = 3/2;
    %Eq. (3.10)
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v) .^ (1 / beta);

    % Eq. (3.9)
    o = step;

end

%lamda 0.5/0.75/1.5
function o = Levy2(d)
    beta = 0.75;
    %Eq. (3.10)
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v) .^ (1 / beta);

    % Eq. (3.9)
    o = step;

end
