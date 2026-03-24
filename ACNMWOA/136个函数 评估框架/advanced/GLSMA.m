% function [Destination_fitness, Convergence_curve]=LASMA(N,Max_FEs,lb,ub,dim,fobj)
function [bestPositions, Convergence_curve] = GLSMA(N, Max_FEs, lb, ub, dim, fobj)

    %初始化参数
    bestPositions = zeros(1, dim); % Destination_position
    Destination_fitness = inf; %对于最大值问题将此改为-inf

    AllFitness = inf * ones(N, 1); %记录所有黏菌的适应度值
    weight = ones(N, dim); %黏菌的适应度值权重

    X = initialization(N, dim, ub, lb); %初始化群体位置
    Convergence_curve = []; %收敛曲线

    it = 1;
    FEs = 1;
    X1 = X;
    z = 0.03; % 参数

    for i = 1:size(X, 1)
        Flag4ub = X(i, :) > ub;
        Flag4lb = X(i, :) < lb;
        X(i, :) = (X(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

        AllFitness(i) = fobj(X(i, :));
        FEs = FEs + 1;
    end

    [SmellOrder, SmellIndex] = sort(AllFitness);
    worstFitness = SmellOrder(size(X, 1)); %初设最后一个是最差的
    bestFitness = SmellOrder(1); %初设第一个是最好的
    S = bestFitness - worstFitness + eps;

    if bestFitness < Destination_fitness
        bestPositions = X1(SmellIndex(1), :);
        Destination_fitness = bestFitness;
    end

    count = 0; %计数
    % 主循环
    while FEs < Max_FEs
        %计算每个黏菌的适应度值
        for i = 1:size(X, 1)

            for j = 2:size(X, 2)

                if i <= (N / 2) %Eq.(2.5)
                    weight(SmellIndex(i), j) = 1 + rand() * log10((bestFitness - SmellOrder(i)) / (S) + 1);
                else
                    weight(SmellIndex(i), j) = 1 - rand() * log10((bestFitness - SmellOrder(i)) / (S) + 1);
                end

            end

        end

        flag = true;

        %更新最佳适应度值和最佳位置
        if bestFitness < Destination_fitness
            bestPositions = X(SmellIndex(1), :);
            Destination_fitness = bestFitness;
            flag = false;
        end

        %%--如果一直没有更新，就使用莱维飞行机制
        if flag
            count = count + 1;
        else
            count = 0;
        end

        if count == 10
            %执行某个搜索
            X_levy = X(i, :) + rand() * sign(rand() -1/2) .* Levy1(dim); % %莱维飞行
            Flag4ub = X_levy > ub;
            Flag4lb = X_levy < lb;
            X(i, :) = X_levy .* (~(Flag4ub + Flag4lb)) + ub .* Flag4ub + lb .* Flag4lb;
        end

        %%----

        a = atanh(- (FEs / Max_FEs) + 1); %Eq.(2.4)
        b = 1 - FEs / Max_FEs;
        XX = zeros(1, dim);
        % 更新每个搜索个体的位置
        for i = 1:N

            if rand < z %Eq.(2.7)
                X(i, :) = (ub - lb) * rand + lb;
            else
                p = tanh(abs(AllFitness(i) - Destination_fitness)); %Eq.(2.2)
                vb = unifrnd(-a, a, 1, dim); %Eq.(2.3)
                vc = unifrnd(-b, b, 1, dim);

                for j = 1:dim
                    r = rand();
                    A = randi([1, N]); % 从种群中随机挑选两个个体
                    B = randi([1, N]);

                    if r < p %Eq.(2.1)
                        X(i, j) = bestPositions(j) + vb(j) * (weight(i, j) * X(A, j) - X(B, j));
                    else
                        X(i, j) = vc(j) * X(i, j);
                    end

                end

            end

            %%%%%%%%%%%%%%%%%%%%%%%高斯变异机制
            X_gaus = X(i, :) * (1 + randn(1)); %randn()是一种产生标准正态分布的随机数或矩阵的函数，

            X_fitness_gaus = fobj(X_gaus);
            X_fitness_s = fobj(X(i, :));
            FEs = FEs + 2;
            X_fitness_comb = [X_fitness_gaus, X_fitness_s];
            [~, m] = min(X_fitness_comb);

            if m == 1
                X(i, :) = X_gaus;
            else
                X(i, :) = X(i, :);
            end

        end

        for i = 1:N
            temp_fit = fobj(X(i, :));
            FEs = FEs + 1;

            if AllFitness(i) > temp_fit
                AllFitness(i) = temp_fit;
                X1(i, :) = X(i, :);
            end

        end

        [SmellOrder, SmellIndex] = sort(AllFitness);
        worstFitness = SmellOrder(size(X, 1));
        bestFitness = SmellOrder(1);
        S = bestFitness - worstFitness + eps;

        if bestFitness < Destination_fitness
            bestPositions = X1(SmellIndex(1), :);
            Destination_fitness = bestFitness;
        end

        Convergence_curve(it) = Destination_fitness;
        it = it + 1;
    end

end

%莱维飞行
function o = Levy1(d)
    beta = 1.5; %beta设置为1.5，与LWOA23相同
    %Eq. (3.10)
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v) .^ (1 / beta);

    % Eq. (3.9)
    o = step;
end
