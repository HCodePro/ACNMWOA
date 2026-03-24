%MFO：飞蛾扑火优化算法
function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = EvolutionWEMFO(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    %输入：
    %N        ：飞蛾数量即行数，影响算法的搜索能力和计算量，一般取20-40，对于较难的问题可以取100
    %MaxFEs：最大评估次数    种群大小*迭代次数=评估次数
    %ub、lb：分别为粒子变量的上下界，如果所有变量的下界都相等，那么可以将lb和ub定义为两个单独的数字。如果变量的上下边界分别不一样，则可以把ub和lb定义为一位数组，例如
    %             lb=[lb1,lb2,...,lbn]，ub=[ub1,ub2,...,ubn]，initialization函数初始化飞蛾种群的时候会判断输入为数字还是数组，从而判断变量的界限是否一致。
    %dim    ：维度，即变量的数量，求解的变量有几个，就设维度为几
    %fobj    ：YourCostFunction 即目标函数，
    %            目标函数是预测值与标签的误差，我们希望找到一个θ使得fobj最小，即误差最小，达到优化的目的。此函数默认为YourCostFunction。可以在单独的文件中定义cost并且加载handle到fobj中

    %输出：
    %Best_flame_pos：优化后的火焰位置，即最优位置
    %Convergence_curve：1行Max_iteration列的矩阵，存放每一次迭代得到的最优结果，对应每一次飞蛾寻找的最优路径，慢慢飞向中心火焰。
    dim = thdim / 2;
    Fbest = -inf;

    Moth_pos = random_initialization(N, dim, ub, lb); %Initialize the positions of moths
    Convergence_curve = zeros(1, MaxFEs); %zeros（a,b）生成全为0的a行b列矩阵
    it = 1; % Iteration 迭代
    Convergence_curve = []; %初始化，空矩阵
    FEs = 0; %评估次数
    %迭代次数和评估次数的区别：种群大小*迭代次数=评估次数。算法会在我们的搜索空间中进行搜索，搜索的点也就是算法中的个体，我们对比的应该是遍历访问了多少区域（点）。
    %所以最后靠谱的应该是“最大评估次数”作为界限
    % Main loop
    s = 0; %qiao

    while FEs < MaxFEs || it <= Iteration %判断是否到达最大评估次数
        Flame_no = round(N - FEs * ((N - 1) / MaxFEs)); % Flame_no 火焰数量 整数  是一个数字，不是数组
        % B = round(A)数组A中每个元素朝最近的方向取整数部分，并返回与A同维的整数数组B  如果A是数字，那返回数字，总之是取整
        %按照等式B=n-g*(n-1)/G	更新火焰数量B， Eq. (3.14) in the paper。
        %火焰数量会随迭代次数增加而逐渐减少，以确保全局和局部搜索之间的平衡。为此，该算法提出了一种自适应调整火焰数量的机制

        for i = 1:size(Moth_pos, 1) % 循环每一个飞蛾 即每一行  i从1一直循环到size(Moth_pos,1)  每次增加1   size(Moth_pos,1）返回Moth_pos的行数
            % Check if moths go out of the search spaceand bring it back，通过变量上下界规范飞蛾个体位置
            Flag4ub = Moth_pos(i, :) > ub; %Moth_pos(i,:)表示Moth_pos的第i行全部数据    一个数组和数字比较大小，结果会产生一个逻辑数组，形式相同，大于则为1，不大于则为0
            Flag4lb = Moth_pos(i, :) < lb; %若飞蛾的位置在最大值和最小值之间，则位置不需要调整，若超出最大值，最回到最大值边界；
            % 若超出最小值，最回到最小值边界
            Moth_pos(i, :) = (Moth_pos(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb; %~表示取反

            if FEs < MaxFEs
                %             FEs=FEs+1;
                % Calculate the fitness of moths
                %             Moth_fitness(1,i)=fobj(Moth_pos(i,:));
                [Moth_fitness(1, i), FEs, Moth_pos(i, :)] = sigle_evaluation(Moth_pos(i, :), dim, thdim, fname, Pxy, FEs);
            else
                break;

            end

        end

        if it == 1 %如果迭代次数为1，Sort the first population of moths，计算飞蛾个体位置的适应度并排序，将排序后对应的个体位置存入火焰矩阵F并将飞蛾个体适应度赋值到火焰适应度矩阵；
            [fitness_sorted I] = sort(Moth_fitness, 'descend'); %sort  对矩阵进行升序排序
            sorted_population = Moth_pos(I, :);

            % Update the flames% 更新火焰
            best_flames = sorted_population;
            best_flame_fitness = fitness_sorted;
            Best_flame_score = fitness_sorted(1); % 更新目前获得的最佳火焰位置

        else %如果不是第一次迭代，则一定已经存在飞蛾种群，下面把当前飞蛾种群和以前的飞蛾种群存到一个容器中，计算适应度值后排序，取最优的前N个，形成新的飞蛾种群。
            double_population = [previous_population; best_flames]; %previous_population 以前的种群   double_population更新的种群
            double_fitness = [previous_fitness best_flame_fitness];

            [double_fitness_sorted I] = sort(double_fitness, 'descend');
            double_sorted_population = double_population(I, :);

            fitness_sorted = double_fitness_sorted(1:N); %计算当前代飞蛾和火焰的适应度并排序，将前n个最优适应度值及其相应位置存入火焰适应度矩阵OF及火焰位置矩阵F作为下一代火焰位置
            sorted_population = double_sorted_population(1:N, :);

            % Update the flames
            best_flames = sorted_population;
            best_flame_fitness = fitness_sorted;
            Best_flame_score = fitness_sorted(1); % 更新目前获得的最佳火焰位置
        end

        %将最优火焰位置存入Fbest

        BB = Best_flame_score;
        Best_flame_pos = sorted_population(1, :);

        previous_population = Moth_pos;
        previous_fitness = Moth_fitness;

        if BB > fitness_sorted(1)
            s = s / 2; %qiao
        end

        s = s + 1; %qiao

        %% 加入双自适应权重，w使算法在前期有较好的全局优化能力，w1使后期有较好的搜索能力
        w1 = (1 - FEs / MaxFEs) ^ (1 - tan(pi * (rand - 0.5)) * s / MaxFEs); %权重w会根据算法陷入局部最优程度呈曲线递减
        w2 = (2 - 2 * FEs / MaxFEs) ^ (1 - tan(pi * (rand - 0.5)) * s / MaxFEs); %权重w1会根据算法陷入局部最优程度呈曲线递减
        a1 = 2 - FEs * ((2) / MaxFEs); % a decreases linearly fron 2 to 0 in Eq. (2.3)
        a = -1 + FEs * ((-1) / MaxFEs);

        % a linearly dicreases(线性分割,划分) from -1 to -2 to calculate t in Eq. (3.12)
        %更新参数b和t,更新飞蛾个体位置

        for i = 1:size(Moth_pos, 1) %对每一只飞蛾循环

            for j = 1:size(Moth_pos, 2) %对列循环

                if i <= Flame_no % Update the position of the moth with respect to its corresponsing flame % 根据相应的火焰更新飞蛾的位置

                    % D in Eq. (3.13) Di = |Fj -  Mi|
                    distance_to_flame = abs(sorted_population(i, j) - Moth_pos(i, j)); %abs返回数组所有元素的绝对值
                    b = 1;
                    t = (a - 1) * rand + 1;

                    % Eq. (3.12) S(Mi,Fj )=Di?e^bt?cos(2πt)+Fj
                    if (FEs / MaxFEs > 0.5)
                        Moth_pos(i, j) = distance_to_flame * exp(b .* t) .* cos(t .* 2 * pi) + w2 * sorted_population(i, j);
                    else
                        Moth_pos(i, j) = w1 * distance_to_flame * exp(b .* t) .* cos(t .* 2 * pi) + sorted_population(i, j);
                    end

                end

                if i > Flame_no % Upudate the position of the moth with respct to one flame %用同一个火焰更新飞蛾的位置

                    % Eq. (3.13)% Eq. (3.13) Di = |Fj -  Mi|
                    distance_to_flame = abs(sorted_population(i, j) - Moth_pos(i, j));
                    b = 1;
                    t = (a - 1) * rand + 1;

                    % Eq. (3.12)    % Eq. (3.12)  S(Mi,Fj )=Di?e^bt?cos(2πt)+Fj
                    if (FEs / MaxFEs > 0.5)
                        Moth_pos(i, j) = w2 * distance_to_flame * exp(b .* t) .* cos(t .* 2 * pi) + sorted_population(Flame_no, j);
                    else
                        Moth_pos(i, j) = distance_to_flame * exp(b .* t) .* cos(t .* 2 * pi) + w1 * sorted_population(Flame_no, j);
                    end

                end

            end

        end

        Convergence_curve(it) = Best_flame_score; %Convergence_curve：1行Max_iteration列的矩阵，存放每一次迭代得到的最优结果，对应每一次飞蛾寻找的最优路径，慢慢飞向中心火焰。

        if Fbest < Best_flame_score
            FE = FEs;
            iter = it;
        end

        Fbest = Best_flame_score;
        Lbest = Best_flame_pos;
        it = it + 1;
    end
