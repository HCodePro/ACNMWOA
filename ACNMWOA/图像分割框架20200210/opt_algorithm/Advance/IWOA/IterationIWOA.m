% The Whale Optimization Algorithm
%% Tubishat, M., et al. (2018). "Improved whale optimization algorithm for feature selection in Arabic sentiment analysis." Applied Intelligence.
% 10.1007/s10489-018-1334-8
%% WOA+OBL+DE
function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = IterationIWOA(fname, SearchAgents_no, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    Fbest = -inf;
    dim = thdim / 2;

    % initialize position vector and score for the leader
    Leader_pos = zeros(1, 2 * dim);
    Leader_score = -inf; %change this to -inf for maximization problems
    FEs = 0;
    Convergence_curve = [];
    t = 1;
    %Initialize the positions of search agents
    Positions = random_initialization(SearchAgents_no, dim, ub, lb);

    %% EOBL
    [P, FEs] = EOBL(Positions, dim, ub, lb, thdim, fname, Pxy, FEs);
    Positions = P;

    for i = 1:size(Positions, 1)
        Flag4ub = Positions(i, :) > ub;
        Flag4lb = Positions(i, :) < lb;
        Positions(i, :) = (Positions(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        [fitness(1, i), FEs, Positions(i, :)] = sigle_evaluation(Positions(i, :), dim, thdim, fname, Pxy, FEs);

        if fitness(1, i) > Leader_score % Change this to > for maximization problem
            Leader_score = fitness(1, i); % Update alpha
            Leader_pos = Positions(i, :);
        end

    end

    % Main loop
    while FEs < MaxFEs || t <= Iteration

        a = 2 - t * ((2) / Iteration); % a decreases linearly fron 2 to 0 in Eq. (2.3)

        % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
        a2 = -1 + t * ((-1) / Iteration);

        % Update the Position of search agents
        for i = 1:size(Positions, 1)
            r1 = rand(); % r1 is a random number in [0,1]
            r2 = rand(); % r2 is a random number in [0,1]

            A = 2 * a * r1 - a; % Eq. (2.3) in the paper
            C = 2 * r2; % Eq. (2.4) in the paper

            b = 1; %  parameters in Eq. (2.5)
            l = (a2 - 1) * rand + 1; %  parameters in Eq. (2.5)

            p = rand(); % p in Eq. (2.6)

            for j = 1:size(Positions, 2)

                if p < 0.5

                    if abs(A) >= 1
                        rand_leader_index = floor(SearchAgents_no * rand() + 1);
                        X_rand = Positions(rand_leader_index, :);
                        D_X_rand = abs(C * X_rand(j) - Positions(i, j)); % Eq. (2.7)
                        Positions(i, j) = X_rand(j) - A * D_X_rand; % Eq. (2.8)

                    elseif abs(A) < 1
                        D_Leader = abs(C * Leader_pos(j) - Positions(i, j)); % Eq. (2.1)
                        Positions(i, j) = Leader_pos(j) - A * D_Leader; % Eq. (2.2)
                    end

                elseif p >= 0.5

                    distance2Leader = abs(Leader_pos(j) - Positions(i, j));
                    % Eq. (2.5)
                    Positions(i, j) = distance2Leader * exp(b .* l) .* cos(l .* 2 * pi) + Leader_pos(j);

                end

            end

        end

        %% Boundary Control
        for i = 1:size(Positions, 1)
            % Return back the search agents that go beyond the boundaries of the search space
            Flag4ub = Positions(i, :) > ub;
            Flag4lb = Positions(i, :) < lb;
            Positions(i, :) = (Positions(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

        end

        %% DE
        for i = 1:size(SearchAgents_no, 1)
            %Mutation
            r = randperm(SearchAgents_no, 3);
            F = rand();
            V(i, :) = Positions(r(1), :) + F * (Positions(r(2), :) - Positions(r(3), :));

            %Cossover
            CR = 0.1;
            jr = randperm(2 * dim, 1);

            for j = 1:2 * dim

                if rand() <= CR || j == jr
                    U(i, j) = V(i, j);
                else
                    U(i, j) = Positions(i, j);
                end

            end

            % Return back the search agents that go beyond the boundaries of the search space
            Flag4ub = U(i, :) > ub;
            Flag4lb = U(i, :) < lb;
            U(i, :) = (U(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            %Selection
            [Ufitness(1, i), FEs, U(i, :)] = sigle_evaluation(U(i, :), dim, thdim, fname, Pxy, FEs);
            [fitness(1, i), FEs, Positions(i, :)] = sigle_evaluation(Positions(i, :), dim, thdim, fname, Pxy, FEs);

            if Ufitness(1, i) > fitness(1, i)
                Positions(i, :) = U(i, :);
                fitness(1, i) = Ufitness(1, i);
            else
                Positions(i, :) = Positions(i, :);
                fitness(1, i) = fitness(1, i);
            end

        end

        [~, Leader_index] = max(fitness);

        if fitness(1, Leader_index) > Leader_score
            Leader_pos = Positions(Leader_index, :);
            Leader_score = fitness(1, Leader_index);
        end

        Convergence_curve(t) = Leader_score;

        if Fbest < Leader_score
            FE = FEs;
            iter = t;
        end

        Fbest = Leader_score;
        Lbest = Leader_pos;
        t = t + 1;

    end
