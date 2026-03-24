% function [BestSol,Convergence_curve]=DE( SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = EvolutionACDE(fname, SearchAgents_no, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    % addpath([pwd,'\statics']);
    dim = thdim / 2;

    Fbest = -inf;
    %% CMA
    lu = [lb; ub];
    Cn = 2 * dim;
    rand('seed', sum(100 * clock));
    CN = Cn;
    maxeval = 300000; % stop criteria
    xmeanw = (lu(1, :) + rand(1, CN) .* (lu(2, :) - lu(1, :)))'; % object parameter start point
    sigma = 0.25; minsigma = 1e-15; maxsigma = max(lu(2, :) - lu(1, :)) / sqrt(Cn); % initial step size, minimal step size
    flginiphase = 1; % Initial phase
    celambda = 4 + floor(3 * log(CN));
    mu = floor(celambda / 2);
    arweights = log((celambda + 1) / 2) - log(1:mu)';
    cc = 4 / (CN + 4); ccov = 2 / (CN + 2 ^ 0.5) ^ 2;
    cs = 4 / (CN + 4); damp = (1 - min(0.7, CN * celambda / maxeval)) / cs + 1;
    B = eye(CN); D = eye(CN); BD = B * D; C = BD * transpose(BD);
    pc = zeros(CN, 1); ps = zeros(CN, 1);
    cw = sum(arweights) / norm(arweights); chiN = CN ^ 0.5 * (1 - 1 / (4 * CN) + 1 / (21 * CN ^ 2));
    counteval = 0; flgstop = 0;
    % Boundary
    Clb = (ones(celambda, 1) * lu(1, :))';
    Cub = (ones(celambda, 1) * lu(2, :))';
    %----------------------------------------------------------------
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
    VarSize = [1, 2 * dim]; % 决策变量矩阵大小 Decision Variables Matrix Size
    beta_min = para(3); % 缩放因子下界 Lower Bound of Scaling Factor
    beta_max = para(4); % 缩放因子上界 Upper Bound of Scaling Factor
    pCR = para(5); %  交叉概率 Crossover Probability
    Lb = ones(1, 2 * dim) .* lb; % 参数取值下界
    Ub = ones(1, 2 * dim) .* ub; % 参数取值上界

    %% 初始化 Initialization
    FEs = 0;
    BestCost = -inf; % 最优值初始化
    BestPos = zeros(1, 2 * dim);
    Position = random_initialization(nPop, dim, ub, lb);
    Cost = zeros(nPop, 1);

    for i = 1:nPop % 遍历每个个体
        [Cost(i), FEs, Position(i, :)] = sigle_evaluation(Position(i, :), dim, thdim, fname, Pxy, FEs);

        if Cost(i) > BestCost % 如果个体目标函数值优于当前最优值
            BestPos = Position(i, :); % 更新最优值
            BestCost = Cost(i);
        end

    end

    Convergence_curve = [];
    it = 1;
    % % Algorithm parameters definition
    Ac = 0.95;
    fr1 = 0.15;
    fr2 = 0.6;
    p1 = 0.2;
    p2 = 0.8;
    tr1 = 0.9;
    tr2 = 0.85;
    tr3 = 0.9;
    %% 主循环 DE Main Loop
    while FEs < MaxIt || it <= Iteration

        newPosition = Position;
        newCost = Cost;

        A_D = squareform(pdist(newPosition, 'squaredeuclidean')); % Eq (4)
        m = tanh(FEs, MaxFEs, [-2, 7]); % Eq (12)
        newPosition_memory = newPosition;
        newCost_memory = newCost;

        for i = 1:nPop
            Dimax = max(A_D(i, :));
            k = floor((1 - FEs / MaxFEs) * nPop) + 1; % Eq (9)
            [~, neighbors] = sort(A_D(i, :));

            % Attraction-Repulsion operator % Eq (6)
            delta_ni = zeros(1, 2 * dim);

            % 假设 epsilon 是您要添加的非常小的数值
            epsilon = 1e-10;

            for j = neighbors(1:k)
                A_I = 1 - (A_D(i, j) / (Dimax + epsilon)); % Eq (7)
                s = sign(newCost(j) - newCost(i)); % Eq (8) % 求max
                delta_ni = delta_ni + Ac * (newPosition_memory(j, :) - newPosition_memory(i, :)) * A_I * s;
            end

            ni = delta_ni / nPop; % Eq (6)
            % Attraction-Repulsion operator % Eq (6)

            % Attraction to best solusion Eq (11)
            if rand < p1
                bi = m * Ac .* (rand(1, 2 * dim) .* BestPos - newPosition_memory(i, :));
            else
                bi = m * Ac .* (BestPos - newPosition_memory(i, :));
            end

            % Attraction to best solusion Eq (10)

            % Local search operators Eq (15)
            if rand < p2

                if rand > 0.5 * FEs / MaxFEs + 0.25
                    u1 = rand(1, 2 * dim) > tr1;
                    ri = u1 .* random('Normal', zeros(1, 2 * dim), fr1 * (1 - FEs / MaxFEs) * (ub - lb)); % Eq (12)
                else
                    u2 = rand(1, 2 * dim) > tr2;
                    w = index_roulette_wheel_selection(newCost, k);
                    Xw = newPosition_memory(w, :);
                    % Eq (14)
                    if rand < 0.5
                        ri = fr2 * u2 .* (1 - FEs / MaxFEs) .* sin(2 * pi * rand(1, 2 * dim)) .* abs(rand(1, 2 * dim) .* Xw - newPosition_memory(i, :));
                    else
                        ri = fr2 * u2 .* (1 - FEs / MaxFEs) .* cos(2 * pi * rand(1, 2 * dim)) .* abs(rand(1, 2 * dim) .* Xw - newPosition_memory(i, :));
                    end

                end

            else
                u3 = rand(1, 2 * dim) > tr3;
                ri = u3 .* (2 * rand(1, 2 * dim) - ones(1, 2 * dim)) .* (ub - lb); % Eq (14)
            end

            % Local search operators Eq (15)

            newPosition(i, :) = newPosition(i, :) + ni + bi + ri; % Eq(16)
        end

        [newPosition, newCost] = memory_operator(newPosition, newCost, newPosition_memory, newCost_memory); % Eq (18)

        for i = 1:nPop
            %Boundary absorption
            Flag4ub = newPosition(i, :) > ub;
            Flag4lb = newPosition(i, :) < lb;
            newPosition(i, :) = (newPosition(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            [newCost(i), FEs, newPosition(i, :)] = sigle_evaluation(newPosition(i, :), dim, thdim, fname, Pxy, FEs);
            %Positive greedy selection mechanism
            if newCost(i) > Cost(i)
                Position(i, :) = newPosition(i, :);
                Cost(i) = newCost(i);

                if newCost(i) > BestCost
                    BestCost = newCost(i);
                    BestPos = newPosition(i, :);
                end

            end

        end

        for i = 1:nPop % 遍历每个个体
            x = newPosition(i, :); % 提取个体位置
            % 随机选择三个个体以备变异使用
            A = randperm(nPop); % 个体顺序重新随机排列
            A(A == i) = []; % 当前个体所排位置腾空（产生变异中间体时当前个体不参与）
            a = A(1);
            b = A(2);
            c = A(3);
            % 变异操作 Mutation
            beta = unifrnd(beta_min, beta_max, VarSize); % 随机产生缩放因子
            y = newPosition(a, :) + beta .* (newPosition(b, :) - newPosition(c, :)); % 产生中间体
            % 防止中间体越界
            y = max(y, Lb);
            y = min(y, Ub);
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

            newPosition(i, :) = z; % 交叉操作之后得到新个体

            [newCost(i), FEs, newPosition(i, :)] = sigle_evaluation(newPosition(i, :), dim, thdim, fname, Pxy, FEs);

            %Positive greedy selection mechanism
            if newCost(i) > Cost(i)
                Position(i, :) = newPosition(i, :);
                Cost(i) = newCost(i);

                if newCost(i) > BestCost
                    BestPos = newPosition(i, :);
                    BestCost = newCost(i);
                end

            end

            %% CMA
            % Generate and evaluate celambda offspring
            arz = randn(CN, celambda);
            arx = xmeanw * ones(1, celambda) + sigma * (BD * arz);

            I = find(arx > Cub);
            arx(I) = 2 * Cub(I) - arx(I);
            aa = find(arx(I) < Clb(I));
            arx(I(aa)) = Clb(I(aa));
            I = find(arx < Clb);
            arx(I) = 2 * Clb(I) - arx(I);
            aa = find(arx(I) > Cub(I));
            arx(I(aa)) = Cub(I(aa));

            U = arx';
            arfitness = zeros(celambda, 1);

            for j = 1:celambda
                Flag4ub = U(j, :) > ub;
                Flag4lb = U(j, :) < lb;
                U(j, :) = (U(j, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
                [fitness, FEs, U(j, :)] = sigle_evaluation(U(j, :), dim, thdim, fname, Pxy, FEs);
                arfitness(j) = fitness;

            end

            counteval = counteval + celambda;
            % Sort by fitness and compute weighted mean in xmeanw
            [arfitness, arindex] = sort(arfitness, 'descend'); % minimization
            xold = xmeanw; % for speed up of Eq. (14)
            xmeanw = arx(:, arindex(1:mu)) * arweights / sum(arweights);
            zmeanw = arz(:, arindex(1:mu)) * arweights / sum(arweights);

            % Adapt covariance matrix
            pc = (1 - cc) * pc + (sqrt(cc * (2 - cc)) * cw / sigma) * (xmeanw - xold); % Eq. (14)

            if ~flginiphase % do not adapt in the initial phase
                C = (1 - ccov) * C + ccov * pc * transpose(pc); % Eq. (15)
            end

            % adapt sigma
            ps = (1 - cs) * ps + (sqrt(cs * (2 - cs)) * cw) * (B * zmeanw); % Eq. (16)
            sigma = sigma * exp((norm(ps) - chiN) / chiN / damp); % Eq. (17)

            % Update B and D from C
            if mod(counteval / celambda, 1 / ccov / CN / 5) < 1
                C = triu(C) + transpose(triu(C, 1)); % enforce symmetry
                [B, D] = eig(C);
                % limit condition of C to 1e14 + 1
                if max(diag(D)) > 1e14 * min(diag(D))
                    tmp = max(diag(D)) / 1e14 - min(diag(D));
                    C = C + tmp * eye(CN); D = D + tmp * eye(CN);
                end

                D = diag(sqrt(diag(D))); % D contains standard deviations now
                BD = B * D; % for speed up only
            end % if mod

            % Adjust minimal step size
            if sigma * min(diag(D)) < minsigma ...
                    | arfitness(1) == arfitness(min(mu + 1, celambda)) ...
                    | xmeanw == xmeanw ...
                    + 0.2 * sigma * BD(:, 1 + floor(mod(counteval / celambda, CN)))
                sigma = 1.4 * sigma;

                % flgstop = 1;
            end

            if sigma > maxsigma
                sigma = maxsigma;
            end

            % Test for end of initial phase
            if flginiphase && counteval / celambda > 2 / cs

                if (norm(ps) - chiN) / chiN < 0.05 % step size is not much too small
                    flginiphase = 0;
                end

            end

            if (arfitness(1) > BestCost)
                BestCost = arfitness(1);
                BestPos = U(arindex(1), :);

            end

        end

        % 保存当前迭代最优个体函数值 Update Best Cost
        Convergence_curve(it) = BestCost;

        if Fbest < BestCost
            FE = FEs;
            iter = it;
        end

        Fbest = BestCost;
        Lbest = BestPos;
        it = it + 1;
    end

end

function [X, F] = memory_operator(X, F, X_memory, F_memory)
    dim = size(X, 2);
    Inx = F_memory > F; % 求max
    Indx = repmat(Inx, 1, dim);
    X = Indx .* X_memory + ~Indx .* X;
    F = Inx .* F_memory + ~Inx .* F;
end

function [y] = tanh(FEs, MaxFEs, range)
    z = 2 * (FEs / MaxFEs * (range(2) - range(1)) + range(1));
    y = 0.5 * ((exp(z) - 1) / (exp(z) + 1) + 1);
end

function [selected_index] = index_roulette_wheel_selection(F, k)
    fitness = F(1:k);
    weights = fitness - min(fitness); % 求max
    weights = cumsum(weights / sum(weights));

    selected_index = roulette_wheel_selection(weights);
end

function [selected_index] = roulette_wheel_selection(weights)
    r = rand();
    selected_index = 1;

    for index = size(weights, 1)

        if r <= weights(index)
            selected_index = index;
            break;
        end

    end

end
