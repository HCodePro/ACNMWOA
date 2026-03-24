function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = IterationSCMINFO(fname, nP, thdim, lb, ub, MaxFEs, Pxy, Iteration)

    %% Initialization
    dim = thdim / 2;
    Fbest = -inf;
    FEs = 0;
    Time = 1;
    Cost = zeros(nP, 1);
    M = zeros(nP, 1);

    X = random_initialization(nP, dim, ub, lb);

    for i = 1:nP
        [fitness, FEs, X(i, :)] = sigle_evaluation(X(i, :), dim, thdim, fname, Pxy, FEs);
        Cost(i) = fitness;
        M(i) = Cost(i);
    end

    [~, ind] = sort(Cost, 'descend');
    Best_X = X(ind(1), :);
    Best_Cost = Cost(ind(1));

    Worst_Cost = Cost(ind(end));
    Worst_X = X(ind(end), :);

    I = randi([2 5]);
    Better_X = X(ind(I), :);
    Better_Cost = Cost(ind(I));

    %% Main Loop of INFO
    while FEs < MaxFEs || Time <= Iteration
        alpha = 2 * exp(-4 * (Time / Iteration)); % Eqs. (5.1) & % Eq. (9.1)

        M_Best = Best_Cost;
        M_Better = Better_Cost;
        M_Worst = Worst_Cost;

        for i = 1:nP

            % Updating rule stage
            del = 2 * rand * alpha - alpha; % Eq. (5)
            sigm = 2 * rand * alpha - alpha; % Eq. (9)

            % Select three random solution
            A1 = randperm(nP);
            A1(A1 == i) = [];
            a = A1(1); b = A1(2); c = A1(3);

            e = 1e-25;
            epsi = e * rand;

            omg = max([M(a) M(b) M(c)]);
            MM = [(M(a) - M(b)) (M(a) - M(c)) (M(b) - M(c))];

            W(1) = cos(MM(1) + pi) * exp(- (MM(1)) / omg); % Eq. (4.2)
            W(2) = cos(MM(2) + pi) * exp(- (MM(2)) / omg); % Eq. (4.3)
            W(3) = cos(MM(3) + pi) * exp(- (MM(3)) / omg); % Eq. (4.4)
            Wt = sum(W);

            WM1 = del .* (W(1) .* (X(a, :) - X(b, :)) + W(2) .* (X(a, :) - X(c, :)) + ... % Eq. (4.1)
                W(3) .* (X(b, :) - X(c, :))) / (Wt + 1) + epsi;

            omg = max([M_Best M_Better M_Worst]);
            MM = [(M_Best - M_Better) (M_Best - M_Better) (M_Better - M_Worst)];

            W(1) = cos(MM(1) + pi) * exp(-MM(1) / omg); % Eq. (4.7)
            W(2) = cos(MM(2) + pi) * exp(-MM(2) / omg); % Eq. (4.8)
            W(3) = cos(MM(3) + pi) * exp(-MM(3) / omg); % Eq. (4.9)
            Wt = sum(W);

            WM2 = del .* (W(1) .* (Best_X - Better_X) + W(2) .* (Best_X - Worst_X) + ... % Eq. (4.6)
                W(3) .* (Better_X - Worst_X)) / (Wt + 1) + epsi;

            % Determine MeanRule
            r = unifrnd(0.1, 0.5);
            MeanRule = r .* WM1 + (1 - r) .* WM2; % Eq. (4)

            if rand < 0.5
                z1 = X(i, :) + sigm .* (rand .* MeanRule) + randn .* (Best_X - X(a, :)) / (M_Best - M(a) + 1);
                z2 = Best_X + sigm .* (rand .* MeanRule) + randn .* (X(a, :) - X(b, :)) / (M(a) - M(b) + 1);
            else % Eq. (8)
                z1 = X(a, :) + sigm .* (rand .* MeanRule) + randn .* (X(b, :) - X(c, :)) / (M(b) - M(c) + 1);
                z2 = Better_X + sigm .* (rand .* MeanRule) + randn .* (X(a, :) - X(b, :)) / (M(a) - M(b) + 1);
            end

            % Vector combining stage
            u = zeros(1, dim * 2);

            for j = 1:dim * 2
                mu = 0.05 * randn;

                if rand < 0.5

                    if rand < 0.5
                        u(j) = z1(j) + mu * abs(z1(j) - z2(j)); % Eq. (10.1)
                    else
                        u(j) = z2(j) + mu * abs(z1(j) - z2(j)); % Eq. (10.2)
                    end

                else
                    u(j) = X(i, j); % Eq. (10.3)
                end

            end

            % Local search stage
            if rand < 0.5
                L = rand < 0.5; v1 = (1 - L) * 2 * (rand) + L; v2 = rand .* L + (1 - L); % Eqs. (11.5) & % Eq. (11.6)
                Xavg = (X(a, :) + X(b, :) + X(c, :)) / 3; % Eq. (11.4)
                phi = rand;
                Xrnd = phi .* (Xavg) + (1 - phi) * (phi .* Better_X + (1 - phi) .* Best_X); % Eq. (11.3)
                Randn = L .* randn(1, dim * 2) + (1 - L) .* randn;

                if rand < 0.5
                    u = Best_X + Randn .* (MeanRule + randn .* (Best_X - X(a, :))); % Eq. (11.1)
                else
                    u = Xrnd + Randn .* (MeanRule + randn .* (v1 * Best_X - v2 * Xrnd)); % Eq. (11.2)
                end

            end

            % Check if new solution go outside the search space and bring them back
            New_X = BC(u, lb, ub);
            [fitness, FEs, New_X] = sigle_evaluation(New_X, dim, thdim, fname, Pxy, FEs);
            New_Cost = fitness;

            if New_Cost > Cost(i)
                X(i, :) = New_X;
                Cost(i) = New_Cost;
                M(i) = Cost(i);

                if Cost(i) > Best_Cost
                    Best_X = X(i, :);
                    Best_Cost = Cost(i);
                end

            end

        end

        % Determine the worst solution
        [~, ind] = sort(Cost, "descend");
        Worst_X = X(ind(end), :);
        Worst_Cost = Cost(ind(end));
        % Determine the better solution
        I = randi([2 5]);
        Better_X = X(ind(I), :);
        Better_Cost = Cost(ind(I));

        % Update Convergence_curve
        Convergence_curve(Time) = Best_Cost;

        if Fbest < Best_Cost
            FE = FEs;
            iter = Time;
        end

        Fbest = Best_Cost;
        Lbest = Best_X;
        Time = Time + 1;

    end

end

function X = BC(X, lb, ub)
    X = fillmissing(X, 'constant', rand);
    Flag4ub = X > ub;
    Flag4lb = X < lb;
    X = (X .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
end
