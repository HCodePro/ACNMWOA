function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = EvolutionWHRIME(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    Fbest = -inf;
    dim = thdim / 2;

    % initialize position
    Best_rime = zeros(1, 2 * dim);
    Best_rime_rate = -inf; %change this to -inf for maximization problems
    Rimepop = random_initialization(N, dim, ub, lb); %Initialize the set of random solutions
    Lb = lb .* ones(1, 2 * dim); % lower boundary
    Ub = ub .* ones(1, 2 * dim); % upper boundary
    % it=1;%Number of iterations
    FEs = 0;
    Time = 1;
    Convergence_curve = [];
    Rime_rates = zeros(1, N); %Initialize the fitness value
    newRime_rates = zeros(1, N);
    W = 5; %Soft-rime parameters, discussed in subsection 4.3.1 of the paper
    %%
    K = 9; %精英数量
    K0 = 8;
    M = zeros(N, 1);
    % flag_K = K;
    %% Calculate the fitness value of the initial position
    for i = 1:N
        %     Rime_rates(1,i)=fobj(Rimepop(i,:));%Calculate the fitness value for each search agent
        [Rime_rates(1, i), FEs, Rimepop(i, :)] = sigle_evaluation(Rimepop(i, :), dim, thdim, fname, Pxy, FEs);
        M(i) = Rime_rates(1, i);
        %     FEs=FEs+1;
        %Make greedy selections
        if Rime_rates(1, i) > Best_rime_rate
            Best_rime_rate = Rime_rates(1, i);
            Best_rime = Rimepop(i, :);
        end

    end

    [~, ind] = sort(Rime_rates, 'descend');
    Best_X = Rimepop(ind(1), :);
    Best_Cost = Rime_rates(ind(1));

    Worst_Cost = Rime_rates(ind(end));
    Worst_X = Rimepop(ind(end), :);

    I = randi([2 5]);
    Better_X = Rimepop(ind(I), :);
    Better_Cost = Rime_rates(ind(I));
    % Main loop
    while FEs < MaxFEs || Time <= Iteration
        %%
        [rate_sort, index] = sort(Rime_rates, 'descend');
        Rime_rates = rate_sort;
        Rimepop = Rimepop(index, :);
        %%
        RimeFactor = (rand - 0.5) * 2 * cos((pi * FEs / (MaxFEs / 10))) * (1 - round(FEs * W / MaxFEs) / W); %Parameters of Eq.(3),(4),(5)
        E = sqrt(FEs / MaxFEs); %Eq.(6)
        newRimepop = Rimepop; %Recording new populations
        normalized_rime_rates = normr(Rime_rates); %Parameters of Eq.(7)
        %% Hierarchical guided strategy
        k = K; %6
        flag = 0;
        flag_K = K;

        while flag < flag_K
            MV = (Rimepop(k, :) + Best_rime) / 2;

            for i = (k + 1):N
                % Determine the better solution
                for j = 1:2 * dim
                    %Soft-rime search strategy
                    r1 = rand();

                    if r1 < E
                        newRimepop(i, j) = Rimepop(k, j) + RimeFactor * ((Ub(j) - Lb(j)) * rand + Lb(j)); %Eq.(3)
                    end

                    %Hard-rime puncture mechanism
                    r2 = rand();

                    if r2 < normalized_rime_rates(i)
                        %                     newRimepop(i,j) = Best_rime(1,j);%Eq.(7)
                        alpha = round(1 + rand());
                        newRimepop(i, j) = Best_rime(1, j) + (Best_rime(1, j) - MV(1, j) * alpha); % 适应度差的个体向优�?��?�?
                    end

                end

                %%
                %             newRime_rates(1,i) = fobj(newRimepop(i,:));
                % Return back the search agents that go beyond the boundaries of the search space
                Flag4ub = newRimepop(i, :) > ub;
                Flag4lb = newRimepop(i, :) < lb;
                newRimepop(i, :) = (newRimepop(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
                [newRime_rates(1, i), FEs, newRimepop(i, :)] = sigle_evaluation(newRimepop(i, :), dim, thdim, fname, Pxy, FEs);
                %             FEs=FEs+1;
                if newRime_rates(1, i) > Rime_rates(1, i)
                    Rime_rates(1, i) = newRime_rates(1, i);
                    M(i) = Rime_rates(1, i);
                    Rimepop(i, :) = newRimepop(i, :);
                    %                 if newRime_rates(1,i) < Best_rime_rate
                    %                    Best_rime_rate=Rime_rates(1,i);
                    %                    Best_rime=Rimepop(i,:);
                    %                 end
                else
                    %% Hierarchical guided strategy
                    M_Best = Best_Cost;
                    M_Better = Better_Cost;
                    M_Worst = Worst_Cost;
                    alpha = 2 * exp(-4 * (FEs / MaxFEs));
                    % Updating rule stage
                    del = 2 * rand * alpha - alpha; % Eq. (5)
                    sigm = 2 * rand * alpha - alpha; % Eq. (9)

                    % Select three random solution
                    A1 = randperm(N);
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

                    WM1 = del .* (W(1) .* (Rimepop(a, :) - Rimepop(b, :)) + W(2) .* (Rimepop(a, :) - Rimepop(c, :)) + ... % Eq. (4.1)
                        W(3) .* (Rimepop(b, :) - Rimepop(c, :))) / (Wt + 1) + epsi;

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
                        z1 = Rimepop(i, :) + sigm .* (rand .* MeanRule) + randn .* (Best_X - Rimepop(a, :)) / (M_Best - M(a) + 1);
                        z2 = Best_X + sigm .* (rand .* MeanRule) + randn .* (Rimepop(a, :) - Rimepop(b, :)) / (M(a) - M(b) + 1);
                    else % Eq. (8)
                        z1 = Rimepop(a, :) + sigm .* (rand .* MeanRule) + randn .* (Rimepop(b, :) - Rimepop(c, :)) / (M(b) - M(c) + 1);
                        z2 = Better_X + sigm .* (rand .* MeanRule) + randn .* (Rimepop(a, :) - Rimepop(b, :)) / (M(a) - M(b) + 1);
                    end

                    % Vector combining stage
                    u = zeros(1, 2 * dim);

                    for j = 1:2 * dim
                        mu = 0.05 * randn;

                        if rand < 0.5

                            if rand < 0.5
                                u(j) = z1(j) + mu * abs(z1(j) - z2(j)); % Eq. (10.1)
                            else
                                u(j) = z2(j) + mu * abs(z1(j) - z2(j)); % Eq. (10.2)
                            end

                        else
                            u(j) = Rimepop(i, j); % Eq. (10.3)
                        end

                    end

                    %                 F_u = fobj(u);
                    % Return back the search agents that go beyond the boundaries of the search space
                    Flag4ub = u > ub;
                    Flag4lb = u < lb;
                    u = (u .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

                    if sum(isnan(u)) == 0
                        [F_u, FEs, u] = sigle_evaluation(u, dim, thdim, fname, Pxy, FEs);
                        %                 FEs=FEs+1;
                        if F_u > Rime_rates(1, i)
                            Rime_rates(1, i) = F_u;
                            M(i) = Rime_rates(1, i);
                            Rimepop(i, :) = u;

                            if F_u > Best_rime_rate
                                Best_rime_rate = F_u;
                                Best_rime = u;
                            end

                        end

                    else
                        break;
                    end

                end

                %             end
            end

            [~, ind] = sort(Rime_rates, 'descend');
            Best_X = Rimepop(ind(1), :);
            Best_Cost = Rime_rates(ind(1));

            Worst_Cost = Rime_rates(ind(end));
            Worst_X = Rimepop(ind(end), :);

            I = randi([2 5]);
            Better_X = Rimepop(ind(I), :);
            Better_Cost = Rime_rates(ind(I));
            %%
            k = k - 1;
            flag = flag + 1;
        end

        K = K0 + floor((1 - FEs / MaxFEs) * (K - K0));
        %%
        Convergence_curve(Time) = Best_rime_rate;

        if Fbest < Best_rime_rate
            FE = FEs;
            iter = Time;
        end

        Fbest = Best_rime_rate;
        Lbest = Best_rime;
        Time = Time + 1;
    end

end
