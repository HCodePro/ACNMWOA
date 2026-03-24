% 单纯�?
% function [Best_universe,Convergence_curve]=IterationMVO(N,MaxFEs,lb,ub,dim,fobj)
function [Fbest, Lbest, FE, MaxFEs, cg_curve, iter] = EvolutionSCMINFO(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    dim = thdim / 2;
    %Two variables for saving the position and inflation rate (fitness) of the best universe
    Best_universe = zeros(1, dim * 2);
    Best_universe_Inflation_rate = -inf;
    Loser_score = inf;
    Loser_pos = zeros(1, dim * 2);
    %Initialize the positions of universes
    % Universes=initialization(N,dim,ub,lb);
    Universes = random_initialization(N, dim, ub, lb);

    %Minimum and maximum of Wormhole Existence Probability (min and max in
    % Eq.(3.3) in the paper
    WEP_Max = 1;
    WEP_Min = 0.2;
    FEs = 0;
    Convergence_curve = [];
    cg_curve = [];

    %Iteration(time) counter
    Time = 1;

    %Main loop
    while FEs < MaxFEs || Time <= Iteration

        %Eq. (3.3) in the paper
        WEP = WEP_Min + FEs * ((WEP_Max - WEP_Min) / MaxFEs);

        %Travelling Distance Rate (Formula): Eq. (3.4) in the paper
        TDR = 1 - ((FEs) ^ (1/6) / (MaxFEs) ^ (1/6));

        %Inflation rates (I) (fitness values)
        Inflation_rates = zeros(1, size(Universes, 1));

        for i = 1:size(Universes, 1)

            %Boundary checking (to bring back the universes inside search
            % space if they go beyoud the boundaries
            Flag4ub = Universes(i, :) > ub;
            Flag4lb = Universes(i, :) < lb;
            Universes(i, :) = (Universes(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            %Calculate the inflation rate (fitness) of universes
            %         Inflation_rates(1,i)=fobj(Universes(i,:));
            if FEs < MaxFEs
                [fitness, FEs, Universes(i, :)] = sigle_evaluation(Universes(i, :), dim, thdim, fname, Pxy, FEs);
            else
                break;
            end

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

        Universes(1, :) = Sorted_universes(1, :);

        %Update the Position of universes
        for i = 2:size(Universes, 1) %Starting from 2 since the firt one is the elite
            Back_hole_index = i;

            for j = 1:size(Universes, 2)
                r1 = rand();

                if r1 < normalized_sorted_Inflation_rates(i)
                    White_hole_index1 = RouletteWheelSelection(-sorted_Inflation_rates); % for maximization problem -sorted_Inflation_rates should be written as sorted_Inflation_rates
                    White_hole_index2 = RouletteWheelSelection(-sorted_Inflation_rates);
                    White_hole_index3 = RouletteWheelSelection(-sorted_Inflation_rates);

                    if White_hole_index1 == -1
                        White_hole_index1 = 1;
                    end

                    if White_hole_index2 == -1
                        White_hole_index2 = 1;
                    end

                    if White_hole_index3 == -1
                        White_hole_index3 = 1;
                    end

                    %Eq. (3.1) in the paper
                    Universes(Back_hole_index, j) = (Sorted_universes(White_hole_index1, j) + Sorted_universes(White_hole_index2, j) + Sorted_universes(White_hole_index3, j)) / 3;
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

            %%%%%%%%%%%%%%%%%%%利用单纯形法更新位置向量%%%%%%%%%%%%%%%%%%%%%%
            %         Pfitness = fobj(Universes(i,:));
            if FEs < MaxFEs
                [fitness, FEs, Universes(i, :)] = sigle_evaluation(Universes(i, :), dim, thdim, fname, Pxy, FEs);
                Pfitness = fitness;
            else
                break;
            end

            if FEs < MaxFEs
                [SimplexPosition, Best_universe_Inflation_rate, Best_universe, FEs] = simplex_han(Best_universe, Best_universe_Inflation_rate, Secondary_pos, Pfitness, Universes(i, :), FEs, dim, thdim, fname, Pxy);
                %         SimplexFitness = fobj(SimplexPosition);

                if FEs < MaxFEs
                    [fitness, FEs, Universes(i, :)] = sigle_evaluation(Universes(i, :), dim, thdim, fname, Pxy, FEs);
                    SimplexFitness = fitness;
                else
                    break;
                end

                if SimplexFitness > Best_universe_Inflation_rate
                    Best_universe_Inflation_rate = SimplexFitness;
                    Best_universe = SimplexPosition;
                end

                if SimplexFitness < Loser_score
                    Loser_score = SimplexFitness;
                end

                if Pfitness > Best_universe_Inflation_rate
                    Best_universe = Universes(i, :);
                    Best_universe_Inflation_rate = Pfitness;
                end

                if Pfitness < Loser_score
                    Loser_score = Pfitness;
                end

                if SimplexFitness > Pfitness
                    Universes(i, :) = SimplexPosition;
                end

                if SimplexFitness < Loser_score
                    Loser_score = SimplexFitness;
                    %Loser_pos = SimplexPosition;
                end

                if SimplexFitness < Best_universe_Inflation_rate && SimplexFitness > Secondary_score
                    Secondary_pos = SimplexPosition;
                    Secondary_score = SimplexFitness;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            else
                break;
            end

        end

        if FEs < MaxFEs

            %Starting Optimizer========================================================
            [bestPoint, bestfitness, allpoint, allfitness] = DE_Process(Universes, Inflation_rates, FEs, ub, lb, Best_universe(1, :), Best_universe_Inflation_rate, dim, thdim, fname, Pxy);
            Best_universe_Inflation_rate = bestfitness;
            Best_universe(1, :) = bestPoint;
        else
            break;
        end

        %Update the convergence curve
        Convergence_curve(Time) = Best_universe_Inflation_rate;

        cg_curve(it) = Best_universe_Inflation_rate;

        if Fbest > Best_universe_Inflation_rate
            FE = FEs;
            iter = it;
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

function [bestPoint, bestfitness, allpoint, allfitness] = DE_Process(allpoint, allfitness, FEs, ub, lb, bestPoint, bestfitness, dim, thdim, fname, Pxy)

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
        z = zeros(size(x)); % 初�?�化一�?新个�?
        j0 = randi([1, numel(x)]); % 产生一�?�?随机数，即选取待交换维度编号？？？

        for j = 1:numel(x) % 遍历每个维度

            if j == j0 || rand <= pCR % 如果当前维度�?待交换维度或者随机�?�率小于交叉概率
                z(j) = y(j); % 新个体当前维度值等于中间体对应维度�?
            else
                z(j) = x(j); % 新个体当前维度值等于当前个体�?�应维度�?
            end

        end

        newpoint = z; % 交叉操作之后得到新个�?
        %     newfitness=fobj(newpoint); % 新个体目标函数�?

        if FEs < MaxFEs
            [fitness, FEs, Universes(i, :)] = sigle_evaluation(Universes(i, :), dim, thdim, fname, Pxy, FEs);
            newfitness = fitness;
        else
            break;
        end

        if newfitness > allfitness(1, i) % 如果新个体优于当前个�?
            allpoint(i, :) = newpoint;
            allfitness(1, i) = newfitness;

            if newfitness > bestfitness
                bestfitness = newfitness;
                bestPoint(1, :) = newpoint;
            end

        end

    end

end

function [Loser_pos, Leader_score, Leader_pos, FEs] = simplex_han(Leader_pos, Leader_score, Secondary_pos, Loser_score, Loser_pos, FEs, dim, thdim, fname, Pxy)
    alpha = 1;
    gamma = 2;
    beta = 0.5;
    Xc = (Leader_pos + Secondary_pos) / 2;
    Xr = Xc + alpha * (Xc - Loser_pos);
    %     XRfitness = fobj(Xr);
    %     FEs = FEs + 1;
    [fitness, FEs, Xr] = sigle_evaluation(Xr, dim, thdim, fname, Pxy, FEs);
    XRfitness = fitness;

    if XRfitness < Leader_score
        Leader_score = XRfitness;
        Leader_pos = Xr;
        Xe = Xc + gamma * (Xr - Xc);
        %         XeFitness = fobj(Xe);
        %         FEs = FEs + 1;
        [fitness, FEs, Xe] = sigle_evaluation(Xe, dim, thdim, fname, Pxy, FEs);
        XeFitness = fitness;

        if XeFitness < Leader_score
            Leader_score = XeFitness;
            Leader_pos = Xe;
            Loser_pos = Xe;
        else
            Loser_pos = Xr;
        end

    end

    if XRfitness > Loser_score
        Xt = Xc + beta * (Loser_pos - Xc);
        %         XtFitness = fobj(Xt);
        %         FEs = FEs + 1;
        [fitness, FEs, Xt] = sigle_evaluation(Xt, dim, thdim, fname, Pxy, FEs);
        XtFitness = fitness;

        if XtFitness < Leader_score
            Leader_score = XtFitness;
            Leader_pos = Xt;
            Loser_pos = Xt;
        else
            Loser_pos = Xr;
        end

    end

    if Leader_score < XRfitness && XRfitness < Loser_score
        Xw = Xc - beta * (Loser_pos - Xc);
        %         XwFitness = fobj(Xw);
        %         FEs = FEs + 1;
        [fitness, FEs, Xw] = sigle_evaluation(Xw, dim, thdim, fname, Pxy, FEs);
        XwFitness = fitness;

        if XwFitness < Leader_score
            Leader_score = XwFitness;
            Leader_pos = Xw;
        end

        if XwFitness < Loser_score
            Loser_pos = Xw;
        else
            Loser_pos = Xr;
        end

    end

end
