% function [BestSol,Convergence_curve]=DE( SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = EvolutionDE(fname, SearchAgents_no, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    % addpath([pwd,'\statics']);
    dim = thdim / 2;
    % fname='renyi2d';
    % SearchAgents_no=50;
    % lb=XVmin;
    % ub=XVmax;
    Fbest = -inf;

    % 参数向量 parameters [n,N_iteration,beta_min,beta_max,pCR]
    % n为种群规模，N_iteration为迭代次数
    % beta_min 缩放因子下界 Lower Bound of Scaling Factor
    % beta_max=0.8; % 缩放因子上界 Upper Bound of Scaling Factor
    % pCR 交叉概率 Crossover Probability
    % 要求输入数据为列向量（矩阵）

    para = [SearchAgents_no, MaxFEs, 0.2, 0.8, 0.2];
    %% 差分进化（DE）算法
    nPop = para(1); % 种群规模 Population Size
    MaxIt = para(2); % 最大迭代次数Maximum Number of Iterations
    nVar = 2 * dim; % 自变量维数，此例需要优化两个参数c和g Number of Decision Variables
    VarSize = [1, 2 * dim]; % 决策变量矩阵大小 Decision Variables Matrix Size
    beta_min = para(3); % 缩放因子下界 Lower Bound of Scaling Factor
    beta_max = para(4); % 缩放因子上界 Upper Bound of Scaling Factor
    pCR = para(5); %  交叉概率 Crossover Probability
    lb = ones(1, 2 * dim) .* lb; % 参数取值下界
    ub = ones(1, 2 * dim) .* ub; % 参数取值上界

    %% 初始化 Initialization
    FEs = 0;
    empty_individual.Position = []; % 种群初始化
    empty_individual.Cost = []; % 种群目标函数值初始化
    % BestSol.Cost=inf; % 最优值初始化
    BestSol.Cost = -inf; % 最优值初始化
    pop = repmat(empty_individual, nPop, 1); % 将保存种群信息的结构体扩展为结构体矩阵，行数等于种群大小

    X = random_initialization(nPop, dim, ub, lb);

    for i = 1:nPop % 遍历每个个体
        %     pop(i).Position=init_individual(lb,ub,2*dim,1);% 随机初始化个体

        pop(i).Position = X(i, :);

        %     pop(i).Cost=fobj(pop(i).Position) ;% 计算个体目标函数值
        %     FEs=FEs+1;

        [pop(i).Cost, FEs, pop(i).Position] = sigle_evaluation(pop(i).Position, dim, thdim, fname, Pxy, FEs);

        if pop(i).Cost > BestSol.Cost % 如果个体目标函数值优于当前最优值
            BestSol = pop(i); % 更新最优值
        end

    end

    BestCost = zeros(MaxIt, 1); % 初始化迭代最优值
    Convergence_curve = [];
    it = 1;
    %% 主循环 DE Main Loop
    while FEs < MaxIt || it <= Iteration

        for i = 1:nPop % 遍历每个个体
            x = pop(i).Position; % 提取个体位置
            % 随机选择三个个体以备变异使用
            A = randperm(nPop); % 个体顺序重新随机排列
            A(A == i) = []; % 当前个体所排位置腾空（产生变异中间体时当前个体不参与）
            a = A(1);
            b = A(2);
            c = A(3);
            % 变异操作 Mutation
            beta = unifrnd(beta_min, beta_max, VarSize); % 随机产生缩放因子
            y = pop(a).Position + beta .* (pop(b).Position - pop(c).Position); % 产生中间体
            % 防止中间体越界
            y = max(y, lb);
            y = min(y, ub);
            % 交叉操作 Crossover
            z = zeros(size(x)); % 初始化一个新个体
            j0 = randi([1, numel(x)]); % 产生一个伪随机数，即选取待交换维度编号？？？

            for j = 1:numel(x) % 遍历每个维度

                if j == j0 || rand <= pCR % 如果当前维度是待交换维度或者随机概率小于交叉概率
                    z(j) = y(j); % 新个体当前维度值等于中间体对应维度值
                else
                    z(j) = x(j); % 新个体当前维度值等于当前个体对应维度值
                end

            end

            NewSol.Position = z; % 交叉操作之后得到新个体
            %         NewSol.Cost=fobj(NewSol.Position); % 新个体目标函数值
            %         FEs=FEs+1;
            if FEs < MaxIt
                [NewSol.Cost, FEs, NewSol.Position] = sigle_evaluation(NewSol.Position, dim, thdim, fname, Pxy, FEs);
            else
                break;
            end

            if NewSol.Cost > pop(i).Cost % 如果新个体优于当前个体
                pop(i) = NewSol; % 更新当前个体

                if pop(i).Cost > BestSol.Cost % 如果当前个体（更新后的）优于最优个体
                    BestSol = pop(i); % 更新最优个体
                end

            end

        end

        % 保存当前迭代最优个体函数值 Update Best Cost
        BestCost(it) = BestSol.Cost;
        Convergence_curve(it) = BestSol.Cost;

        %     display(['Iteration ' num2str(it) ': Best fitness = ' num2str(Convergence_curve(it))]);

        if Fbest < BestSol.Cost
            FE = FEs;
            iter = it;
        end

        Fbest = BestSol.Cost;
        Lbest = BestSol.Position;
        it = it + 1;
    end

    % Positionbest=BestSol.Position;
    % bestCVaccuarcy=BestSol.Cost;
    % Fbest=bestCVaccuarcy;
    % Lbest=Positionbest;

end
