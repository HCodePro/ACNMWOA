% The Whale Optimization Algorithm
%% Chaos+CSO+NMs
function [Leader_pos, Convergence_curve] = ACNMWOA(SearchAgents_no, MaxFEs, lb, ub, dim, fobj)

    % initialize position vector and score for the leader
    Leader_pos = zeros(1, dim);
    Leader_score = inf; %change this to -inf for maximization problems

    %Initialize the positions of search agents
    Positions = initialization(SearchAgents_no, dim, ub, lb);
    Lb = ones(1, dim) .* lb;
    Ub = ones(1, dim) .* ub;
    Convergence_curve = [];
    fitness = zeros(SearchAgents_no, 1);
    FEs = 0;
    t = 1;

    % Main loop
    while FEs < MaxFEs

        a = 2 - FEs * ((2) / MaxFEs); % a decreases linearly fron 2 to 0 in Eq. (2.3)

        % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
        a2 = -1 + FEs * ((-1) / MaxFEs);
        %% start chaos

        chFitness = zeros(size(fitness));
        chPositions = zeros(size(Positions));
        %% start chaos
        a_1 = 1;
        index = 4;

        for i = 1:SearchAgents_no
            x = rand();
            O = chaos(index, dim, a_1, x, lb, ub);

            chPositions(i, :) = Positions(i, :) .* O;
            Flag4ub = chPositions(i, :) > Ub;
            Flag4lb = chPositions(i, :) < Lb;
            chPositions(i, :) = (chPositions(i, :) .* (~(Flag4ub + Flag4lb))) + Ub .* Flag4ub + Lb .* Flag4lb;
            chFitness(i) = fobj(chPositions(i, :));
            FEs = FEs + 1;

            if chFitness(i) < fitness(i)
                Positions(i, :) = chPositions(i, :);
                fitness(i) = chFitness(i);

                if chFitness(i) < Leader_score
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
                fitness(i) = fobj(Positions(i, :));
                FEs = FEs + 1;
                % Update the leader
                if fitness(i) < Leader_score % Change this to > for maximization problem
                    Leader_score = fitness(i); % Update alpha
                    Leader_pos = Positions(i, :);
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
            [Positions, fitness, FEs] = CC(Positions, fitness, dim, Lb, Ub, fobj, FEs);
            % i_1 = floor(SearchAgents_no * rand() + 1);
            i_1 = 27;

            for j = 1:size(Positions, 2)

                if p < 0.5

                    if abs(A) >= 1
                        % 寻找猎物
                        rand_leader_index = floor(SearchAgents_no * rand() + 1);
                        X_rand = Positions(rand_leader_index, :);
                        D_X_rand = abs(C * X_rand(j) - Positions(i_1, j)); % Eq. (2.7)
                        Positions(i_1, j) = X_rand(j) - A * D_X_rand; % Eq. (2.8)

                    elseif abs(A) < 1
                        % 收缩环绕更新位置 气泡攻击
                        D_Leader = abs(C * Leader_pos(j) - Positions(i_1, j)); % Eq. (2.1)
                        Positions(i_1, j) = Leader_pos(j) - A * D_Leader; % Eq. (2.2)
                    end

                elseif p >= 0.5

                    % Positions(i_1, j) = s * Levy(1) * (Positions(i_1, j) - Leader_pos(j)) + Positions(i_1, j);
                    % 螺旋式攻击更新位置 气泡攻击
                    distance2Leader = abs(Leader_pos(j) - Positions(i_1, j));
                    % Eq. (2.5)
                    Positions(i_1, j) = distance2Leader * exp(b .* l) .* cos(l .* 2 * pi) + Leader_pos(j);

                end

            end

        end

        options = optimset('MaxFunEvals', floor(MaxFEs * 0.2));
        [x, fval, ~, output] = fminsearchbnd(fobj, Leader_pos, Lb, Ub, options);

        if fval < Leader_score
            Leader_score = fval;
            Leader_pos = x;
        end

        FEs = FEs + output.funcCount;
        Convergence_curve(t) = Leader_score;
        t = t + 1;

    end

end
