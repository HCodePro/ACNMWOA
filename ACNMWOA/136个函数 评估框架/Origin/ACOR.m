function [best_solution, convergence_curve] = ACOR(SearchAgents_no, MaxFEs, lb, ub, dim, fobj)
    FEs = 0;
    it = 1;
    rho = 0.5; % 信息素挥发速率
    num_ants = SearchAgents_no;
    num_iterations = MaxFEs / SearchAgents_no;
    pheromone_initial = 0.1; % 初�?�信�?素�?

    % 初�?�化
    best_solution = zeros(1, dim);
    best_fitness = inf;
    pheromone = pheromone_initial * ones(1, dim);
    convergence_curve = [];

    for iter = 1:num_iterations
        solutions = lb + (ub - lb) * rand(num_ants, dim); % 初�?�化蚂蚁的解决方�?

        % 计算每只蚂蚁的适应度�?
        fitness_values = zeros(num_ants, 1);

        for ant = 1:num_ants
            fitness_values(ant) = fobj(solutions(ant, :));
            FEs = FEs + 1;
        end

        % 更新最佳解
        [min_fitness, min_index] = min(fitness_values);

        if min_fitness < best_fitness
            best_fitness = min_fitness;
            best_solution = solutions(min_index, :);
        end

        % 更新信息�?
        for ant = 1:num_ants
            delta_pheromone = 1 / fitness_values(ant); % 更新信息素的增量
            pheromone = (1 - rho) * pheromone + delta_pheromone;
        end

        % 更新收敛曲线
        convergence_curve(it) = best_fitness;
        it = it + 1;
    end

end
