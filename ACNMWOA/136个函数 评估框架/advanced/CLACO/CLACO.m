function [Leader_pos, Convergence_curve] = CLACO(SearchAgents_no, MaxFEs, lb, ub, dim, fobj)
    %% Problem Definition
    solution = [1 dim]; %初始化一个解决方法的大小
    %% 初始化相关参数
    FEs = 0;
    it = 1;
    k = 10; % 档案袋大小k
    m = SearchAgents_no; % 新生成个数
    q = 0.5;
    ibslo = 1;
    %% 初始化
    % 初始化一个个体
    empty_individual.Position = [];
    empty_individual.fitness = [];
    % 初始化种群
    pop = repmat(empty_individual, k, 1);

    for i = 1:k
        pop(i).Position = unifrnd(lb, ub, solution);
        Flag4ub = pop(i).Position > ub; %返回大于ub的逻辑值
        Flag4lb = pop(i).Position < lb;
        pop(i).Position = (pop(i).Position .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        pop(i).fitness = fobj(pop(i).Position);
        FEs = FEs + 1;
    end

    % 按照函数值排序
    [~, SortOrder] = sort([pop.fitness]);
    pop = pop(SortOrder);
    Bestsolution = pop(1); % 保存一次迭代的最佳值的解决方案
    % Convergence_curve=zeros(MaxFEs,1);            % 保留每一次迭代的最佳值
    w = 1 / (sqrt(2 * pi) * q * k) * exp(-0.5 * (((1:k) - 1) / (q * k)) .^ 2); %计算权重
    p = w / sum(w); % 计算概率
    k1 = 0;

    wMax = 0.9;
    wMin = 0.4;
    epsilon = 30;
    lambda = 0.01;
    a = 2;
    r = zeros(1, dim);
    rMax = zeros(1, dim);
    theta = zeros(1, dim);
    beta = 1.5;
    deltaU = ((gamma(1 + beta) * sin(pi * beta / 2)) / (gamma((1 + beta) / 2) * beta * (2 ^ ((beta - 1) / 2)))) ^ (1 / beta);
    deltaV = 1;
    gs_best = zeros(1, dim);

    %% 主循环
    while FEs <= MaxFEs
        % 原k个方案
        s = zeros(k, dim);

        for l = 1:k
            s(l, :) = pop(l).Position(1, :);
        end

        % 新m个方案
        newpop = repmat(empty_individual, m, 1);

        for t = 1:m
            % 轮盘赌选择高斯函数
            for j = 1:k
                p1 = cumsum(p);

                if p1(j) >= rand
                    l = j;
                    break;
                end

            end

            % 初始化新的解决方案
            newpop(t).Position = zeros(solution);
            % 生成新的解决方案
            for i = 1:dim
                % 标准差计算
                D = 0;

                for r = 1:k
                    D = D + abs(s(l, i) - s(r, i));
                end

                sigma = ibslo * D / (k - 1);
                % 根据高斯核函数产生高斯随机值
                %             newpop(t).Position(i)=normrnd(s(l,i),sigma);
                newpop(t).Position(i) = s(l, i) + sigma * randn;
            end

            Flag4ub = newpop(t).Position > ub; %返回大于ub的逻辑值
            Flag4lb = newpop(t).Position < lb;
            newpop(t).Position = (newpop(t).Position .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            % Evaluation
            newpop(t).fitness = fobj(newpop(t).Position);
            FEs = FEs + 1;

            %% 柯西变异
            %         if(tan(pi*(rand-0.5))<(1-FEs/MaxFEs))
            if FEs / MaxFEs > 0.5
                r = cauchyrnd(0, 1, dim);
                t1 = randperm(dim);
                cauchy = r(t1(1), :);
                X_meaution(1, :) = s(l, :) + rand(1, dim) .* cauchy;
                Flag4ub = X_meaution > ub; %返回大于ub的逻辑值
                Flag4lb = X_meaution < lb;
                X_meaution = (X_meaution .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
                % Evaluation
                X_meaution_fitness = fobj(X_meaution);
                FEs = FEs + 1;

                if (X_meaution_fitness < pop(1).fitness)
                    pop(1).fitness = X_meaution_fitness;
                    pop(1).Position = X_meaution;
                end

            end

        end

        % 将新生成的m个解加入k个原始方案中
        pop = [pop; newpop];
        % 将所有解排序
        [~, SortOrder] = sort([pop.fitness]);
        pop = pop(SortOrder);
        % 删除m个解
        pop = pop(1:k);
        % 更新最佳值的解决方案
        Bestsolution = pop(1);

        for i = 1:k
            X(i, :) = pop(i).Position;
        end

        %% optimal individual based on greedy levy variation
        for j = 1:dim
            gs_best = Bestsolution.Position;
            r(j) = abs(Bestsolution.Position(j) - (1 / k) * sum(X(:, j)));
            rMax(j) = max(X(:, j)) - min(X(:, j));
            theta(j) = exp((-epsilon * FEs / MaxFEs) * (1 - r(j) / rMax(j)));
            levy = normrnd(0, deltaU ^ 2) / (abs(normrnd(0, deltaV ^ 2)) ^ (1 / beta));
            gs_best(j) = gs_best(j) + theta(j) * levy * Bestsolution.Position(j);
            gs_best_fit = fobj(gs_best);
            FEs = FEs + 1;

            if gs_best_fit < Bestsolution.fitness
                Bestsolution.Position(j) = gs_best(j);
                Bestsolution.fitness = gs_best_fit;
                pop(1).Position(j) = gs_best(j);
                pop(1).fitness = gs_best_fit;
            end

        end

        % 保存每次迭代的最佳函数值
        Convergence_curve(it) = Bestsolution.fitness;
        Leader_pos = Bestsolution.Position;
        % 显示每一次的迭代结果信息
        %     display(['Iteration ' num2str(it) ': Best fitness = ' num2str(Convergence_curve(it))]);
        k1 = 1;
        it = it + 1;
    end

end
