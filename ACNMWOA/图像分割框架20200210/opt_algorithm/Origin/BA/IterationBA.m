% Main programs starts here
% function [best,Convergence_curve]=BA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = IterationBA(fname, SearchAgents_no, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    dim = thdim / 2;
    Fbest = -inf;
    n = SearchAgents_no;
    A = rand(1, n) + ones(1, n); % Loudness  (constant or decreasing)
    r = rand(1, n); % Pulse rate (constant or decreasing)
    % This frequency range determines the scalings
    % You should change these values if necessary
    Qmin = 0; % Frequency minimum
    Qmax = 2; % Frequency maximum

    % 2*dimension of the search variables
    d = 2 * dim; % Number of 2*dimensions
    % Lower limit/bounds/ a vector
    Lb = lb .* ones(1, d);
    % Upper limit/bounds/ a vector
    Ub = ub .* ones(1, d);
    % Initializing arrays
    Q = zeros(n, 1); % Frequency
    v = zeros(n, d); % Velocities
    FEs = 0;
    % Initialize the population/solutions
    Sol = random_initialization(n, dim, ub, lb);

    for i = 1:n
        %   Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
        %   Fitness(i)=fobj(Sol(i,:));
        %   FEs=FEs+1;
        [Fitness(i), FEs, Sol(i, :)] = sigle_evaluation(Sol(i, :), dim, thdim, fname, Pxy, FEs);

    end

    % Find the initial best solution
    [fmin, I] = min(Fitness);
    best = Sol(I, :);

    Convergence_curve = [];

    t = 0;

    % Main loop
    while FEs < MaxFEs || t <= Iteration
        % for t=1:N_gen,
        % Loop over all bats/solutions
        for i = 1:n
            Q = Qmin + (Qmin - Qmax) * rand;
            v(i, :) = v(i, :) + (Sol(i, :) - best) * Q;
            S(i, :) = Sol(i, :) + v(i, :);
            % Apply simple bounds/limits
            Sol(i, :) = simplebounds(Sol(i, :), Lb, Ub);
            % Pulse rate
            if rand > r(i)
                % The factor 0.001 limits the step sizes of random walks
                S(i, :) = best + A(i) * randn(1, d);
            end

            % Evaluate new solutions
            S(i, :) = simplebounds(S(i, :), Lb, Ub);
            %                 FEs=FEs+1;
            %                 Fnew=fobj(S(i,:));
            [Fnew, FEs, S(i, :)] = sigle_evaluation(S(i, :), dim, thdim, fname, Pxy, FEs);

            % Update if the solution improves, or not too loud
            if (Fnew >= Fitness(i))
                Sol(i, :) = S(i, :);
                Fitness(i) = Fnew;
            end

            if (Fnew > fmin) && (rand < A(i))
                A(i) = 0.9 * A(i);
                r(i) = r(i) * (1 - exp((-0.9) * t));
            end

            % Update the current best solution
            if Fnew >= fmin
                best = S(i, :);
                fmin = Fnew;
            end

        end

        t = t + 1;
        Convergence_curve(t) = fmin;

        if Fbest < fmin
            FE = FEs;
            iter = t;
        end

        Fbest = fmin;
        Lbest = best;
        %         display(['Iteration ' num2str(t) ': Best fitness = ' num2str(Convergence_curve(t))]);

    end

end

% Application of simple limits/bounds
function s = simplebounds(s, Lb, Ub)
    % Apply the lower bound vector
    ns_tmp = s;
    I = ns_tmp < Lb;
    ns_tmp(I) = Lb(I);

    % Apply the upper bound vector
    J = ns_tmp > Ub;
    ns_tmp(J) = Ub(J);
    % Update this new move
    s = ns_tmp;
end
