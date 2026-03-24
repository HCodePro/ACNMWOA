% CCDE
%% Clustering Center-based Differential Evolution 10.1109/CEC55065.2022.9870429

function [bestPositions, Convergence_curve] = CCDE(SearchAgents_no, MaxFEs, lb, ub, dim, fobj)
    % 参数向量 parameters [n,N_iteration,beta_min,beta_max,pCR]
    % n为种群规模，N_iteration为迭代次数
    % beta_min 缩放因子下界 Lower Bound of Scaling Factor
    % beta_max=0.8; % 缩放因子上界 Upper Bound of Scaling Factor
    % pCR 交叉概率 Crossover Probability
    % 要求输入数据为列向量（矩阵）
    para = [SearchAgents_no, MaxFEs, 0.2, 0.8, 0.9];
    %% 差分进化（DE）算法
    nPop = para(1); % 种群规模 Population Size
    MaxIt = para(2); % 最大迭代次数Maximum Number of Iterations
    nVar = dim; % 自变量维数，此例需要优化两个参数c和g Number of Decision Variables
    VarSize = [1, dim]; % 决策变量矩阵大小 Decision Variables Matrix Size
    beta_min = para(3); % 缩放因子下界 Lower Bound of Scaling Factor
    beta_max = para(4); % 缩放因子上界 Upper Bound of Scaling Factor
    pCR = para(5); %  交叉概率 Crossover Probability
    lb = ones(1, dim) .* lb; % 参数取值下界
    ub = ones(1, dim) .* ub; % 参数取值上界
    %% 初始化 Initialization
    FEs = 0;
    N = para(1);
    bestPositions = zeros(1, dim);
    Destination_fitness = inf;
    AllFitness = inf * ones(N, 1); %record the fitness of all positions
    X = initialization(N, dim, ub, lb);
    Nc = 5; % 簇的个数
    SC = (N - Nc) / Nc;

    for i = 1:nPop % 遍历每个个体
        %     pop(i).Position=init_individual(lb,ub,dim,1);% 随机初始化个体
        %     pop(i).Cost=fobj(pop(i).Position) ;% 计算个体目标函数值
        AllFitness(i) = fobj(X(i, :));
        FEs = FEs + 1;

        if AllFitness(i) < Destination_fitness % 如果个体目标函数值优于当前最优值
            bestPositions = X(i, :); % 更新最优值
            Destination_fitness = AllFitness(i);
        end

    end

    BestCost = zeros(MaxIt, 1); % 初始化迭代最优值
    Convergence_curve = [];
    it = 1;
    %% 主循环 DE Main Loop
    while FEs < MaxIt

        for i = 1:(nPop - Nc) % 遍历除簇个数以外的每个个体
            x = X(i, :); % 提取个体位置
            % 随机选择三个个体以备变异使用
            A = randperm(nPop); % 个体顺序重新随机排列
            A(A == i) = []; % 当前个体所排位置腾空（产生变异中间体时当前个体不参与）
            a = A(1);
            b = A(2);
            c = A(3);
            % 变异操作 Mutation
            beta = unifrnd(beta_min, beta_max, VarSize); % 随机产生缩放因子
            y = X(a, :) + beta .* (X(b, :) - X(c, :)); % 产生中间体
            % 防止中间体越界
            y = max(y, lb);
            y = min(y, ub);
            % 交叉操作 Crossover
            z = zeros(size(x)); % 初始化一个新个体
            j0 = randi([1, numel(x)]); % 产生一个伪随机数，即选取待交换维度编号

            for j = 1:numel(x) % 遍历每个维度

                if j == j0 || rand <= pCR % 如果当前维度是待交换维度或者随机概率小于交叉概率
                    z(j) = y(j); % 新个体当前维度值等于中间体对应维度值
                else
                    z(j) = x(j); % 新个体当前维度值等于当前个体对应维度值
                end

            end

            Fz = fobj(z); % 新个体目标函数值
            FEs = FEs + 1;

            if Fz < AllFitness(i) % 如果新个体优于当前个体
                X(i, :) = z; % 更新当前个体
                AllFitness(i) = Fz;

                if Fz < Destination_fitness % 如果当前个体（更新后的）优于最优个体
                    bestPositions = z; % 更新最优个体
                    Destination_fitness = Fz;
                end

            end

        end

        %% cluster
        [fitvalue, fitorder] = sort(AllFitness);
        X = X(fitorder, :);
        AllFitness = fitvalue;

        for t = 1:Nc
            start_value = (t - 1) * SC + 1;
            end_value = t * SC;
            X_c = 0;

            for p = start_value:end_value
                X_c = X_c + X(p, :);
            end

            X_c = X_c / SC;
            X(t + (N - Nc), :) = X_c;
            AllFitness(t + (N - Nc)) = fobj(X_c);
            FEs = FEs + 1;
        end

        [minvalue, minindex] = min(AllFitness);

        if minvalue < Destination_fitness
            Destination_fitness = minvalue;
            bestPositions = X(minindex, :);
        end

        % 保存当前迭代最优个体函数值 Update Best Cost
        BestCost(it) = Destination_fitness;
        Convergence_curve(it) = Destination_fitness;
        it = it + 1;
    end

end

function V = ControlBoundD(V, L, U) % 防止越界

    for i = 1:size(V, 1)
        Flag4ub = V(i, :) > U;
        Flag4lb = V(i, :) < L;
        V(i, :) = (V(i, :) .* (~(Flag4ub + Flag4lb))) + U .* Flag4ub + L .* Flag4lb;
    end

end
