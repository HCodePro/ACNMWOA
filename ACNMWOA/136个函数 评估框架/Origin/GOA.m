% The Grasshopper Optimization Algorithm
function [Convergence_curve] = GOA(N, MaxFEs, lb, ub, dim, fobj)

    flag = 0;

    if size(ub, 1) == 1
        ub = ones(1, dim) * ub;
        lb = ones(1, dim) * lb;
    end

    if (rem(dim, 2) ~= 0) % this algorithm should be run with a even number of variables. This line is to handle odd number of variables
        dim = dim + 1;
        ub = [ub, 100];
        lb = [lb, -100];
        flag = 1;
    end

    FEs = 0;
    %Initialize the population of grasshoppers
    GrassHopperPositions = initialization(N, dim, ub, lb);
    GrassHopperFitness = zeros(1, N);
    ub = ub';
    lb = lb';

    fitness_history = zeros(N, MaxFEs);
    position_history = zeros(N, MaxFEs, dim);
    Convergence_curve = [];
    Trajectories = zeros(N, MaxFEs);

    cMax = 1;
    cMin = 0.00004;
    %Calculate the fitness of initial grasshoppers

    for i = 1:size(GrassHopperPositions, 1)

        if flag == 1
            GrassHopperFitness(1, i) = fobj(GrassHopperPositions(i, 1:end - 1));
            FEs = FEs + 1;
        else
            GrassHopperFitness(1, i) = fobj(GrassHopperPositions(i, :));
            FEs = FEs + 1;
        end

        fitness_history(i, 1) = GrassHopperFitness(1, i);
        position_history(i, 1, :) = GrassHopperPositions(i, :);
        Trajectories(:, 1) = GrassHopperPositions(:, 1);
    end

    [sorted_fitness, sorted_indexes] = sort(GrassHopperFitness);

    % Find the best grasshopper (target) in the first population
    for newindex = 1:N
        Sorted_grasshopper(newindex, :) = GrassHopperPositions(sorted_indexes(newindex), :);
    end

    TargetPosition = Sorted_grasshopper(1, :);
    TargetFitness = sorted_fitness(1);

    % Main loop
    l = 2; % Start from the second iteration since the first iteration was dedicated to calculating the fitness of antlions

    while FEs < MaxFEs

        c = cMax - FEs * ((cMax - cMin) / MaxFEs); % Eq. (2.8) in the paper

        for i = 1:size(GrassHopperPositions, 1)
            temp = GrassHopperPositions';

            for k = 1:2:dim
                S_i = zeros(2, 1);

                for j = 1:N

                    if i ~= j
                        Dist = distance(temp(k:k + 1, j), temp(k:k + 1, i)); % Calculate the distance between two grasshoppers

                        r_ij_vec = (temp(k:k + 1, j) - temp(k:k + 1, i)) / (Dist + eps); % xj-xi/dij in Eq. (2.7)
                        xj_xi = 2 + rem(Dist, 2); % |xjd - xid| in Eq. (2.7)

                        s_ij = ((ub(k:k + 1) - lb(k:k + 1)) * c / 2) * S_func(xj_xi) .* r_ij_vec; % The first part inside the big bracket in Eq. (2.7)
                        S_i = S_i + s_ij;
                    end

                end

                S_i_total(k:k + 1, :) = S_i;

            end

            X_new = c * S_i_total' + (TargetPosition); % Eq. (2.7) in the paper
            GrassHopperPositions_temp(i, :) = X_new';
        end

        % GrassHopperPositions
        GrassHopperPositions = GrassHopperPositions_temp;

        for i = 1:size(GrassHopperPositions, 1)
            % Relocate grasshoppers that go outside the search space
            Tp = GrassHopperPositions(i, :) > ub';
            Tm = GrassHopperPositions(i, :) < lb';
            GrassHopperPositions(i, :) = (GrassHopperPositions(i, :) .* (~(Tp + Tm))) + ub' .* Tp + lb' .* Tm;

            % Calculating the objective values for all grasshoppers
            if flag == 1
                GrassHopperFitness(1, i) = fobj(GrassHopperPositions(i, 1:end - 1));
                FEs = FEs + 1;
            else
                GrassHopperFitness(1, i) = fobj(GrassHopperPositions(i, :));
                FEs = FEs + 1;
            end

            fitness_history(i, l) = GrassHopperFitness(1, i);
            position_history(i, l, :) = GrassHopperPositions(i, :);

            Trajectories(:, l) = GrassHopperPositions(:, 1);

            % Update the target
            if GrassHopperFitness(1, i) < TargetFitness
                TargetPosition = GrassHopperPositions(i, :);
                TargetFitness = GrassHopperFitness(1, i);
            end

        end

        Convergence_curve(l) = TargetFitness;

        l = l + 1;
    end

    Convergence_curve(1) = Convergence_curve(2);

    if (flag == 1)
        TargetPosition = TargetPosition(1:dim - 1);
    end

end

function Positions = initialization(SearchAgents_no, dim, ub, lb)

    Boundary_no = size(ub, 2); % numnber of boundaries
    % If the boundaries of all variables are equal and user enter
    % a signle number for both ub and lb
    if Boundary_no == 1
        Positions = (ub - lb) .* rand(SearchAgents_no, dim) + lb * ones(SearchAgents_no, dim);
    end

    % If each variable has a different lb and ub
    if Boundary_no > 1

        for i = 1:dim
            ub_i = ub(i);
            lb_i = lb(i);
            Positions(:, i) = rand(SearchAgents_no, 1) .* (ub_i - lb_i) + lb_i;
        end

    end

end

function d = distance(a, b)
    d = sqrt((a(1) - b(1)) ^ 2 + (a(2) - b(2)) ^ 2);
end

function o = S_func(r)
    f = 0.5;
    l = 1.5;
    o = f * exp(-r / l) - exp(-r); % Eq. (2.3) in the paper
end
