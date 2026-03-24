% Particle Swarm Optimization
% function [bestPos,cg_curve]=PSO(N,MaxFEs,lb,ub,2*dim,fobj)
function [Fbest, Lbest, FE, MaxFEs, cg_curve, iter] = IterationPSO(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)
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
    % pos=initialization(noP,dim,ub,lb);
    pos = random_initialization(noP, dim, ub, lb);

    for i = 1:noP
        pBestScore(i) = -inf;
    end

    % Initialize gBestScore for a minimization problem
    gBestScore = -inf;

    it = 1;
    FEs = 0;

    % Main loop
    while FEs < MaxFEs || it <= Iteration

        % Return back the particles that go beyond the boundaries of the search
        % space
        %      Flag4ub=pos(i,:)>ub;
        %      Flag4lb=pos(i,:)<lb;
        %      pos(i,:)=(pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        for i = 1:size(pos, 1)
            %Calculate objective function for each particle

            Flag4ub = pos(i, :) > ub;
            Flag4lb = pos(i, :) < lb;
            pos(i, :) = (pos(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            %             fitness=fobj(pos(i,:));
            %             FEs=FEs+1;
            [fitness, FEs, pos(i, :)] = sigle_evaluation(pos(i, :), dim, thdim, fname, Pxy, FEs);

            if (pBestScore(i) < fitness)
                pBestScore(i) = fitness;
                pBest(i, :) = pos(i, :);
            end

            if (gBestScore < fitness)
                gBestScore = fitness;
                gBest = pos(i, :);
            end

        end

        %Update the W of PSO
        %     w=wMax-l*((wMax-wMin)/iter);
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

        cg_curve(it) = gBestScore;

        if Fbest < gBestScore
            FE = FEs;
            iter = it;
        end

        Fbest = gBestScore;
        bestPos = gBest;
        Lbest = bestPos;
        %     display(['Iteration ' num2str(it) ': Best fitness = ' num2str(cg_curve(it))]);
        it = it + 1;
    end

end
