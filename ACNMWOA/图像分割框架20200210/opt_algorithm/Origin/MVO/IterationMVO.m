function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = IterationMVO(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)

    %Two variables for saving the position and inflation rate (fitness) of the best universe
    dim = thdim / 2;
    Best_universe = zeros(1, 2 * dim);
    Fbest = -inf;
    Best_universe_Inflation_rate = -inf;
    Loser_score = inf;
    Loser_pos = zeros(1, dim);

    %Initialize the positions of universes
    Universes = random_initialization(N, dim, ub, lb);

    %Minimum and maximum of Wormhole Existence Probability (min and max in
    % Eq.(3.3) in the paper
    WEP_Max = 1;
    WEP_Min = 0.2;
    FEs = 0;
    Convergence_curve = [];

    %Iteration(time) counter
    Time = 1;

    %Main loop
    while FEs < MaxFEs || Time <= Iteration

        %Eq. (3.3) in the paper
        WEP = WEP_Min + Time * ((WEP_Max - WEP_Min) / Iteration);

        %Travelling Distance Rate (Formula): Eq. (3.4) in the paper
        TDR = 1 - ((Time) ^ (1/6) / (Iteration) ^ (1/6));

        %Inflation rates (I) (fitness values)
        Inflation_rates = zeros(1, size(Universes, 1));

        for i = 1:size(Universes, 1)

            %Boundary checking (to bring back the universes inside search
            % space if they go beyoud the boundaries
            Flag4ub = Universes(i, :) > ub;
            Flag4lb = Universes(i, :) < lb;
            Universes(i, :) = (Universes(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            %Calculate the inflation rate (fitness) of universes
            % Inflation_rates(1,i)=fobj(Universes(i,:));
            % FEs=FEs+1;
            [fitness, FEs, Universes(i, :)] = sigle_evaluation(Universes(i, :), dim, thdim, fname, Pxy, FEs);
            Inflation_rates(1, i) = fitness;
            %Elitism
            if Inflation_rates(1, i) > Best_universe_Inflation_rate
                Best_universe_Inflation_rate = Inflation_rates(1, i);
                Best_universe = Universes(i, :);
            end

            if Inflation_rates(1, i) < Loser_score
                Loser_score = Inflation_rates(1, i);
                Loser_pos = Universes(i, :);
            end

        end

        [sortFitness, index] = sort(Inflation_rates, 'descend');
        Secondary_score = sortFitness(N - 1);
        Secondary_pos = Universes(index(N - 1), :);
        [sorted_Inflation_rates, sorted_indexes] = sort(Inflation_rates);

        for newindex = 1:N
            Sorted_universes(newindex, :) = Universes(sorted_indexes(newindex), :);
        end

        %Normaized inflation rates (NI in Eq. (3.1) in the paper)
        normalized_sorted_Inflation_rates = normr(sorted_Inflation_rates);

        Universes(1, :) = Sorted_universes(N, :);

        %Update the Position of universes
        for i = 1:size(Universes, 1) %Starting from 2 since the firt one is the elite
            Back_hole_index = i;

            for j = 1:size(Universes, 2)
                r1 = rand();

                if r1 < normalized_sorted_Inflation_rates(i)
                    White_hole_index = RouletteWheelSelection(-sorted_Inflation_rates); % for maximization problem -sorted_Inflation_rates should be written as sorted_Inflation_rates

                    if White_hole_index == -1
                        White_hole_index = 1;
                    end

                    %Eq. (3.1) in the paper

                    Universes(Back_hole_index, j) = Sorted_universes(White_hole_index, j);
                end

                if (size(lb, 2) == 1)
                    %Eq. (3.2) in the paper if the boundaries are all the same
                    r2 = rand();

                    if r2 < WEP
                        r3 = rand();

                        if r3 < 0.5
                            Universes(i, j) = Best_universe(1, j) + TDR * ((ub - lb) * rand + lb);
                        end

                        if r3 > 0.5
                            Universes(i, j) = Best_universe(1, j) - TDR * ((ub - lb) * rand + lb);
                        end

                    end

                end

                if (size(lb, 2) ~= 1)
                    %Eq. (3.2) in the paper if the upper and lower bounds are
                    %different for each variables
                    r2 = rand();

                    if r2 < WEP
                        r3 = rand();

                        if r3 < 0.5
                            Universes(i, j) = Best_universe(1, j) + TDR * ((ub(j) - lb(j)) * rand + lb(j));
                        end

                        if r3 > 0.5
                            Universes(i, j) = Best_universe(1, j) - TDR * ((ub(j) - lb(j)) * rand + lb(j));
                        end

                    end

                end

            end

        end

        %Update the convergence curve

        Convergence_curve(Time) = Best_universe_Inflation_rate;

        if Fbest < Best_universe_Inflation_rate
            FE = FEs;
            iter = Time;
        end

        Fbest = Best_universe_Inflation_rate;
        Lbest = Best_universe;
        %Print the best universe details after every 50 iterations
        %     if mod(Time,50)==0
        %         display(['At iteration ', num2str(Time), ' the best universes fitness is ', num2str(Best_universe_Inflation_rate)]);
        %     end
        Time = Time + 1;
    end

end

function [bestPoint, bestfitness, allpoint, allfitness, FEs] = DE_Process(allpoint, allfitness, FEs, ub, lb, bestPoint, bestfitness, dim, thdim, fname, Pxy)

    beta_min = 0.2; % 缩放因子下界 Lower Bound of Scaling Factor
    beta_max = 0.8; % 缩放因子上界 Upper Bound of Scaling Factor
    pCR = 0.2; %  交叉概率 Crossover Probability
    VarSize = [1, size(allpoint, 2)];

    for i = 1:size(allpoint, 1)
        x = allpoint(i, :); % 提取�?体位�?
        % 随机选择三个�?体以备变异使�?
        A = randperm(size(allpoint, 1)); % �?体顺序重新随机排�?
        A(A == i) = []; % 当前�?体所排位�?腾空（产生变异中间体时当前个体不参与�?
        a = A(1);
        b = A(2);
        c = A(3);
        % 变异操作 Mutation
        beta = unifrnd(beta_min, beta_max, VarSize); % 随机产生缩放因子
        y = allpoint(a, :) + beta .* (allpoint(b, :) - allpoint(c, :)); % 产生�?间体
        % 防�??�?间体越界
        y = max(y, lb);
        y = min(y, ub);

        % 交叉操作 Crossover
        z = zeros(1, size(x, 2)); % 初�?�化一�?新个�?
        j0 = randi([1, numel(x)]); % 产生一�?�?随机数，即选取待交换维度编号？？？

        for j = 1:numel(x) % 遍历每个维度

            if j == j0 || rand <= pCR % 如果当前维度�?待交换维度或者随机�?�率小于交叉概率
                z(j) = y(j); % 新个体当前维度值等于中间体对应维度�?
            else
                z(j) = x(j); % 新个体当前维度值等于当前个体�?�应维度�?
            end

        end

        newpoint(1, :) = z; % 交叉操作之后得到新个�?
        %     newfitness=fobj(newpoint); % 新个体目标函数�?
        % FEs=FEs+1;
        Flag4ub = newpoint(1, :) > ub;
        Flag4lb = newpoint(1, :) < lb;
        newpoint(1, :) = (newpoint(1, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        [fitness, FEs, newpoint(i, :)] = sigle_evaluation(newpoint(1, :), dim, thdim, fname, Pxy, FEs);
        FEs = FEs;

        if fitness > bestfitness
            allpoint(i, :) = newpoint(1, :);
            allfitness(1, i) = fitness;
            bestfitness = fitness;
            bestPoint(1, :) = newpoint(1, :);
        end

    end

end

%%%%%%%%%%%%%%%%%%%单纯形法%%%%%%%%%%%%%%%%%%%%
function [Loser_pos, Leader_score, Leader_pos, FEs] = simplex_han(lb, ub, dim, thdim, Leader_pos, Leader_score, Secondary_pos, Loser_score, Loser_pos, FEs, fname, Pxy)
    alpha = 1;
    gamma = 2;
    beta = 0.5;
    Xc = (Leader_pos + Secondary_pos) / 2;
    Xr = Xc + alpha * (Xc - Loser_pos);
    Flag4ub = Xr > ub;
    Flag4lb = Xr < lb;
    Xr = (Xr .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
    [XRfitness, FEs, Xr] = sigle_evaluation(Xr, dim, thdim, fname, Pxy, FEs);
    % XRfitness = fobj(Xr);
    % FEs = FEs + 1;
    if XRfitness > Leader_score
        Leader_score = XRfitness;
        Leader_pos = Xr;
        Xe = Xc + gamma * (Xr - Xc);
        Flag4ub = Xe > ub;
        Flag4lb = Xe < lb;
        Xe = (Xe .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        [XeFitness, FEs, Xe] = sigle_evaluation(Xe, dim, thdim, fname, Pxy, FEs);
        %     XeFitness = fobj(Xe);
        %     FEs = FEs + 1;
        if XeFitness > Leader_score
            Leader_score = XeFitness;
            Leader_pos = Xe;
            Loser_pos = Xe;
        else
            Loser_pos = Xr;
        end

    end

    if XRfitness < Loser_score
        Xt = Xc + beta * (Loser_pos - Xc);
        Flag4ub = Xt > ub;
        Flag4lb = Xt < lb;
        Xt = (Xt .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        [XtFitness, FEs, Xt] = sigle_evaluation(Xt, dim, thdim, fname, Pxy, FEs);
        %     XtFitness = fobj(Xt);
        %     FEs = FEs + 1;
        if XtFitness > Leader_score
            Leader_score = XtFitness;
            Leader_pos = Xt;
            Loser_pos = Xt;
        else
            Loser_pos = Xr;
        end

    end

    if Leader_score > XRfitness && XRfitness > Loser_score
        Xw = Xc - beta * (Loser_pos - Xc);
        Flag4ub = Xw > ub;
        Flag4lb = Xw < lb;
        Xw = (Xw .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        [XwFitness, FEs, Xw] = sigle_evaluation(Xw, dim, thdim, fname, Pxy, FEs);
        %     XwFitness = fobj(Xw);
        %     FEs = FEs + 1;
        if XwFitness > Leader_score
            Leader_score = XwFitness;
            Leader_pos = Xw;
        end

        if XwFitness > Loser_score
            Loser_pos = Xw;
        else
            Loser_pos = Xr;
        end

    end

end
