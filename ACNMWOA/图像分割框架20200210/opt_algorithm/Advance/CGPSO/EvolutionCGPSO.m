% Particle Swarm Optimization
function [Fbest, Lbest, FE, MaxFEs, cg_curve, iter] = EvolutionCGPSO(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    %% Tsung-Ying Sun, Chan-Cheng Liu, Shang-Jeng Tsai, Sheng-Ta Hsieh, and Kan-Yuan Li; Cluster Guide Particle Swarm Optimization(CGPSO) for Underdetermined Blind Source Separation with Advanced Conditions;
    %% IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL. 15, NO. 6, DECEMBER 2011
    dim = thdim / 2;
    %PSO Infotmation
    Fbest = -inf;
    Vmax = 6;
    noP = N;
    wMax = 0.9;
    wMin = 0.2;
    c1 = 2;
    c2 = 2;

    % Initializations
    vel = zeros(noP, 2 * dim);
    pBestScore = zeros(noP);
    pBest = zeros(noP, 2 * dim);
    gBest = zeros(1, 2 * dim);
    cg_curve = [];

    % Random initialization for agents.
    pos = random_initialization(noP, dim, ub, lb);

    for i = 1:noP
        pBestScore(i) = -inf;
    end

    % Initialize gBestScore for a minimization problem
    gBestScore = -inf;

    FEs = 0;
    it = 1;

    while FEs < MaxFEs || it <= Iteration

        % Return back the particles that go beyond the boundaries of the search
        % space

        for i = 1:size(pos, 1)
            %Calculate objective function for each particle
            Flag4ub = pos(i, :) > ub;
            Flag4lb = pos(i, :) < lb;
            pos(i, :) = (pos(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            if FEs < MaxFEs
                [fitness, FEs, pos(i, :)] = sigle_evaluation(pos(i, :), dim, thdim, fname, Pxy, FEs);

                if (pBestScore(i) < fitness)
                    pBestScore(i) = fitness;
                    pBest(i, :) = pos(i, :);
                end

                if (gBestScore < fitness)
                    gBestScore = fitness;
                    gBest = pos(i, :);
                end

            else
                break;
            end

        end

        w = 1;
        %Update the Velocity and Position of particles
        for i = 1:size(pos, 1)

            for j = 1:size(pos, 2)
                vel(i, j) = w * vel(i, j) + c1 * rand() * (pBest(i, j) - pos(i, j)) + c2 * rand() * (gBest(j) - pos(i, j));

                if (vel(i, j) > Vmax)
                    vel(i, j) = Vmax;
                end

                if (vel(i, j) <- Vmax)
                    vel(i, j) = -Vmax;
                end

                pos(i, j) = pos(i, j) + vel(i, j);
            end

        end

        %%%=================CLS=================================
        if FEs / MaxFEs < 0.8
            %         setCan = (MaxFEs-FEs+1)/MaxFEs;
            setCan = 1 - power(abs((FEs - 1) / FEs), 1500);

            x = rand();

            if (x ~= 0.25 && x ~= 0.5 && x ~= 0.75)
                ch(1) = x;
            end

            for ant = 1:(noP)
                ch(ant + 1) = 4 * ch(ant) * (1 - ch(ant));
                CH(ant, :) = lb + ch(ant) * (ub - lb); %ub´ó
                V = (1 - setCan) * gBest + setCan * CH(ant);
                Flag4ub = V > ub;
                Flag4lb = V < lb;
                V = (V .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

                if FEs < MaxFEs
                    [FitnessV, FEs, V] = sigle_evaluation(V, dim, thdim, fname, Pxy, FEs);

                    if (FitnessV > gBestScore)
                        gBestScore = FitnessV;
                        gBest = V;
                        break;
                    end

                else
                    break;
                end

            end

        else
            setCan = (MaxFEs - FEs + 1) / MaxFEs;
            %
            ch = randn(1, noP);

            for ant = 1:(noP)
                CH(ant, :) = lb + ch(ant) * (ub - lb); %ub´ó
                V = (1 - setCan) * gBest + setCan * CH(ant);
                Flag4ub = V > ub;
                Flag4lb = V < lb;
                V = (V .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

                if FEs < MaxFEs
                    [FitnessV, FEs, V] = sigle_evaluation(V, dim, thdim, fname, Pxy, FEs);

                    if (FitnessV > gBestScore)
                        gBestScore = FitnessV;
                        gBest = V;
                        break;
                    end

                else
                    break;
                end

            end

        end

        %%%=======================================================
        cg_curve(it) = gBestScore;
        bestPos = gBest;

        if Fbest < gBestScore
            FE = FEs;
            iter = it;
        end

        Fbest = gBestScore;
        Lbest = gBest;
        it = it + 1;
    end

end
