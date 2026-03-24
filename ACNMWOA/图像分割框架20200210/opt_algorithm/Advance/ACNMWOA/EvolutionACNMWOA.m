% The Whale Optimization Algorithm
% function [Leader_pos,Convergence_curve]=WOA(SearchAgents_no,MaxFEs,lb,ub,2*dim,fobj)
function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = EvolutionACNMWOA(fname, SearchAgents_no, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    disp('CCNMWOA is now tackling your problem')
    dim = thdim / 2;
    Fbest = -inf;
    % initialize position vector and score for the leader
    % Leader_pos=zeros(1,dim);
    Leader_pos = zeros(1, 2 * dim);
    Leader_score = -inf; %change this to -inf for maximization problems

    %Initialize the positions of search agents
    % Positions=initialization(SearchAgents_no,dim,ub,lb);
    Positions = random_initialization(SearchAgents_no, dim, ub, lb);
    fitness = zeros(SearchAgents_no, 1);
    Lb = ones(1, 2 * dim) .* lb;
    Ub = ones(1, 2 * dim) .* ub;
    Convergence_curve = [];
    FEs = 0;
    t = 1;

    % Main loop
    while FEs < MaxFEs || t <= Iteration

        a = 2 - FEs * ((2) / MaxFEs); % a decreases linearly fron 2 to 0 in Eq. (2.3)

        % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
        a2 = -1 + FEs * ((-1) / MaxFEs);
        % %% start chaos
        % % Ratio of the left iteration with the Max iteration.
        % setCan = (MaxFEs - FEs + 1) / MaxFEs;
        % % To assign a random value to x that don't equal to 0.25, 0.5, 0.75
        % x = rand();

        % while (~(x ~= 0.25 && x ~= 0.5 && x ~= 0.75))
        %     x = rand();
        % end

        % ch(1) = x;

        % for ant = 1:size(Positions, 1)
        %     ch(ant + 1) = 4 * ch(ant) * (1 - ch(ant));
        %     CH(ant, :) = lb + ch(ant) * (ub - lb); %ub大
        %     V = (1 - setCan) * Leader_pos + setCan * CH(ant);
        %     Flag4ub = V > ub;
        %     Flag4lb = V < lb;
        %     V = (V .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        %     [FitnessV, FEs, V] = sigle_evaluation(V, dim, thdim, fname, Pxy, FEs);
        %     % FitnessV = fobj(V); %计算适应度值
        %     % FEs = FEs + 1;

        %     if (FitnessV > Leader_score)
        %         Leader_score = FitnessV;
        %         Leader_pos = V;
        %         break;
        %     end

        % end

        % %% end chaos
        %% start chaos
        chFitness = zeros(size(fitness));
        chPositions = zeros(size(Positions));
        %% 设置混沌映射的种类
        a_1 = 1;
        index = 4;

        for i = 1:SearchAgents_no
            x = rand();
            O = chaos(index, 2 * dim, a_1, x, lb, ub);

            chPositions(i, :) = Positions(i, :) .* O;
            % chFitness(i) = fobj(chPositions(i, :));
            % FEs = FEs + 1;
            Flag4ub = chPositions(i, :) > Ub;
            Flag4lb = chPositions(i, :) < Lb;
            chPositions(i, :) = (chPositions(i, :) .* (~(Flag4ub + Flag4lb))) + Ub .* Flag4ub + Lb .* Flag4lb;
            [chFitness(i), FEs, chPositions(i, :)] = sigle_evaluation(chPositions(i, :), dim, thdim, fname, Pxy, FEs);

            if chFitness(i) > fitness(i)
                Positions(i, :) = chPositions(i, :);
                fitness(i) = chFitness(i);

                if chFitness(i) > Leader_score
                    Leader_score = fitness(i);
                    Leader_pos = Positions(i, :);
                end

            end

        end

        %% end chaos
        % Update the Position of search agents
        for i_1 = 1:size(Positions, 1)

            for i = 1:size(Positions, 1)

                % Return back the search agents that go beyond the boundaries of the search space
                Flag4ub = Positions(i, :) > Ub;
                Flag4lb = Positions(i, :) < Lb;
                Positions(i, :) = (Positions(i, :) .* (~(Flag4ub + Flag4lb))) + Ub .* Flag4ub + Lb .* Flag4lb;

                % Calculate objective function for each search agent
                %         fitness=fobj(Positions(i,:));
                %         FEs=FEs+1;
                if FEs < MaxFEs
                    [fitness(i), FEs, Positions(i, :)] = sigle_evaluation(Positions(i, :), dim, thdim, fname, Pxy, FEs);

                    % Update the leader
                    if fitness(i) > Leader_score % Change this to > for maximization problem
                        Leader_score = fitness(i); % Update alpha
                        Leader_pos = Positions(i, :);
                    end

                else
                    break;
                end

            end

            r1 = rand(); % r1 is a random number in [0,1]
            r2 = rand(); % r2 is a random number in [0,1]

            A = 2 * a * r1 - a; % Eq. (2.3) in the paper
            C = 2 * r2; % Eq. (2.4) in the paper

            b = 1; %  parameters in Eq. (2.5)
            l = (a2 - 1) * rand + 1; %  parameters in Eq. (2.5)

            p = rand(); % p in Eq. (2.6)
            %% 纵横交叉策略CSO
            [Positions, fitness, FEs] = CC(Positions, fitness, 2 * dim, Lb, Ub, FEs, dim, thdim, fname, Pxy);
            % i_1 = floor(SearchAgents_no * rand() + 1);
            i_1 = 27;

            for j = 1:2 * dim

                if p < 0.5

                    if abs(A) >= 1
                        rand_leader_index = floor(SearchAgents_no * rand() + 1);
                        X_rand = Positions(rand_leader_index, :);
                        D_X_rand = abs(C * X_rand(j) - Positions(i_1, j)); % Eq. (2.7)
                        Positions(i_1, j) = X_rand(j) - A * D_X_rand; % Eq. (2.8)

                    elseif abs(A) < 1
                        D_Leader = abs(C * Leader_pos(j) - Positions(i_1, j)); % Eq. (2.1)
                        Positions(i_1, j) = Leader_pos(j) - A * D_Leader; % Eq. (2.2)
                    end

                elseif p >= 0.5

                    distance2Leader = abs(Leader_pos(j) - Positions(i_1, j));
                    % Eq. (2.5)
                    Positions(i_1, j) = distance2Leader * exp(b .* l) .* cos(l .* 2 * pi) + Leader_pos(j);

                end

            end

        end

        %% 如何对单纯形的方法进行修改
        options = optimset('MaxFunEvals', floor(MaxFEs * 0.2));
        [x, fval, ~, output, params] = fminsearchbnd(Leader_pos, Lb, Ub, options, dim, thdim, fname, Pxy, FEs);
        FEs = params.FEs;

        if fval > Leader_score
            Leader_score = fval;
            Leader_pos = x;
        end

        %FEs = FEs + output.funcCount;
        Convergence_curve(t) = Leader_score;
        %     display(['Iteration ' num2str(t) ': Best fitness = ' num2str(Convergence_curve(t))]);
        if Fbest < Leader_score
            FE = FEs;
            iter = t;
        end

        Fbest = Leader_score;
        Lbest = Leader_pos;
        t = t + 1;
    end
