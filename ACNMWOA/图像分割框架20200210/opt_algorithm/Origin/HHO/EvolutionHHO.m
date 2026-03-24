function [Fbest, Lbest, FE, MaxFES, CNVG, iter] = EvolutionHHO(fname, N, thdim, lb, ub, MaxFES, Pxy, Iteration)
    dim = thdim / 2;
    Fbest = -inf;

    % initialize the location and Energy of the rabbit
    Rabbit_Location = zeros(1, thdim);
    Rabbit_Energy = -inf;

    %Initialize the locations of Harris' hawks
    X = random_initialization(N, dim, ub, lb);

    CNVG = [];
    fes = 0; % FES counter

    t = 1; % Loop counter

    while fes < MaxFES || t <= Iteration

        for i = 1:size(X, 1)
            X(i, :) = BC(X(i, :), lb, ub);
            [newfitness, fes, X(i, :)] = sigle_evaluation(X(i, :), dim, thdim, fname, Pxy, fes);
            fitness = newfitness;

            % Update the location of Rabbit
            if fitness > Rabbit_Energy
                Rabbit_Energy = fitness;
                Rabbit_Location = X(i, :);
            end

        end

        c = 2 * (1 - (fes / MaxFES)); % factor to show the decreaing energy of rabbit
        % Update the location of Harris' hawks
        for i = 1:size(X, 1)

            Escaping_Energy = c * (2 * rand() - 1); % escaping energy of rabbit

            if abs(Escaping_Energy) >= 1
                %% Exploration:
                % Harris' hawks search for the rabbit, observing the area to find a rabbit
                % Harris' hawks perch randomly based on 2 strategy:

                r = rand();
                rand_Hawk_index = floor(N * rand() + 1);
                X_rand = X(rand_Hawk_index, :);

                if r < 0.5
                    % perch based on other family members
                    X(i, :) = X_rand - rand() * abs(X_rand - 2 * rand() * X(i, :));
                elseif r > 0.5
                    % perch on a random tall tree (random site inside group's home range)
                    X(i, :) = (Rabbit_Location(1, :) - mean(X)) - rand() * ((ub - lb) * rand + lb);
                end

            elseif abs(Escaping_Energy) < 1
                %% Exploitation:
                % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit

                %% phase 1: surprise pounce (seven kills): hawks go to kill when rabbit is surprised
                % surprise pounce (seven kills): multiple, short rapid dives by different hawks

                r = rand(); % probablity of each event % r>0.5 : rabbit cannot escape -- r<0.5 : rabbit can escape

                if r > 0.5 && abs(Escaping_Energy) < 0.5 % Hard besiege % When rabbit cannot or do not try to escape
                    X(i, :) = (Rabbit_Location) - Escaping_Energy * abs(Rabbit_Location - X(i, :));
                end

                if r > 0.5 && abs(Escaping_Energy) > 0.5 % Soft besiege % When rabbit cannot or do not try to escape
                    Jump_strength = 2 * (1 - rand()); % random jump strength of the rabbit
                    X(i, :) = (Rabbit_Location - X(i, :)) - Escaping_Energy * abs(Jump_strength * Rabbit_Location - X(i, :));
                end

                %% phase 2: performing team rapid dives (leapfrog movements) when prey escapes
                if r < 0.5 && abs(Escaping_Energy) > 0.5, % Soft besiege % rabbit try to escape by many zigzag deceptive motions

                    Jump_strength = 2 * (1 - rand());
                    X1 = Rabbit_Location - Escaping_Energy * abs(Jump_strength * Rabbit_Location - X(i, :));

                    X1 = BC(X1, lb, ub);
                    X(i, :) = BC(X(i, :), lb, ub);
                    [newfitness1, fes, X1] = sigle_evaluation(X1, dim, thdim, fname, Pxy, fes);
                    [newfitness2, fes, X(i, :)] = sigle_evaluation(X(i, :), dim, thdim, fname, Pxy, fes);

                    if newfitness1 > newfitness2 % improved move?
                        X(i, :) = X1;
                    else % hawks perform levy-based short rapid dives around the rabbit
                        X2 = Rabbit_Location - Escaping_Energy * abs(Jump_strength * Rabbit_Location - X(i, :)) + rand(1, thdim) .* Levy(thdim);
                        X2 = BC(X2, lb, ub);
                        [newfitness3, fes, X2] = sigle_evaluation(X2, dim, thdim, fname, Pxy, fes);

                        if (newfitness3 > newfitness2), % improved move?
                            X(i, :) = X2;
                        end

                    end

                end

                if r < 0.5 && abs(Escaping_Energy) < 0.5, % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                    % hawks try to decrease their average location with the rabbit
                    Jump_strength = 2 * (1 - rand());
                    X1 = Rabbit_Location - Escaping_Energy * abs(Jump_strength * Rabbit_Location - mean(X));

                    X1 = BC(X1, lb, ub);
                    X(i, :) = BC(X(i, :), lb, ub);
                    [newfitness1, fes, X1] = sigle_evaluation(X1, dim, thdim, fname, Pxy, fes);
                    [newfitness2, fes, X(i, :)] = sigle_evaluation(X(i, :), dim, thdim, fname, Pxy, fes);

                    if newfitness1 > newfitness2 % improved move?
                        X(i, :) = X1;
                    else % Perform levy-based short rapid dives around the rabbit
                        X2 = Rabbit_Location - Escaping_Energy * abs(Jump_strength * Rabbit_Location - mean(X)) + rand(1, thdim) .* Levy(thdim);
                        X2 = BC(X2, lb, ub);
                        [newfitness3, fes, X2] = sigle_evaluation(X2, dim, thdim, fname, Pxy, fes);

                        if (newfitness3 > newfitness2), % improved move?
                            X(i, :) = X2;
                        end

                    end

                end

                %%
            end

        end

        CNVG(t) = Rabbit_Energy;

        if Fbest < Rabbit_Energy
            FE = fes;
            iter = t;
        end

        Fbest = Rabbit_Energy;
        Lbest = Rabbit_Location;
        t = t + 1;

    end

end

function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);
    u = randn(1, d) * sigma; v = randn(1, d); step = u ./ abs(v) .^ (1 / beta);
    o = step;
end

function X = BC(X, lb, ub)
    X = fillmissing(X, 'constant', rand);
    Flag4ub = X > ub;
    Flag4lb = X < lb;
    X = (X .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
end
