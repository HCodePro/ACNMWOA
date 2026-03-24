function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = IterationRIME(fname, SearchAgents_no, thdim, lb, ub, MaxFEs, Pxy, Iteration)
    % disp('RIME is now tackling your problem')
    dim = thdim / 2;

    Fbest = -inf; %å…¨éƒ¨è¿?ä»£æ??ä¼˜å??
    % initialize position
    Best_rime = zeros(1, 2 * dim);
    %Best_rime_rate=inf;%change this to -inf for maximization problems
    Best_rime_rate = -inf; %ä¸?æ¬¡è¿­ä»£æœ€ä¼˜å??
    Rimepop = random_initialization(SearchAgents_no, dim, ub, lb); %Initialize the set of random solutions
    Lb = lb .* ones(1, 2 * dim); % lower boundary
    Ub = ub .* ones(1, 2 * dim); % upper boundary
    % it=1;%Number of iterations
    FEs = 0;
    Time = 1;
    Convergence_curve = [];

    Rime_rates = zeros(1, SearchAgents_no); %Initialize the fitness value
    newRime_rates = zeros(1, SearchAgents_no);
    W = 5; %Soft-rime parameters, discussed in subsection 4.3.1 of the paper
    %Calculate the fitness value of the initial position
    for i = 1:SearchAgents_no
        %Rime_rates(1,i)=fobj(Rimepop(i,:));%Calculate the fitness value for each search agent
        [Rime_rates(1, i), FEs, Rimepop(i, :)] = sigle_evaluation(Rimepop(i, :), dim, thdim, fname, Pxy, FEs);
        %FEs=FEs+1;
        %Make greedy selections
        %if Rime_rates(1,i)<Best_rime_rate
        if Rime_rates(1, i) > Best_rime_rate
            Best_rime_rate = Rime_rates(1, i);
            Best_rime = Rimepop(i, :);
        end

    end

    % Main loop
    while FEs < MaxFEs || Time <= Iteration

        RimeFactor = (rand - 0.5) * 2 * cos((pi * Time / (Iteration / 10))) * (1 - round(Time * W / Iteration) / W); %Parameters of Eq.(3),(4),(5)
        E = sqrt(Time / Iteration); %Eq.(6)
        newRimepop = Rimepop; %Recording new populations
        normalized_rime_rates = normr(Rime_rates); %Parameters of Eq.(7)

        for i = 1:SearchAgents_no

            for j = 1:2 * dim
                %Soft-rime search strategy
                r1 = rand();

                if r1 < E
                    newRimepop(i, j) = Best_rime(1, j) + RimeFactor * ((Ub(j) - Lb(j)) * rand + Lb(j)); %Eq.(3)

                end

                %Hard-rime puncture mechanism
                r2 = rand();

                if r2 < normalized_rime_rates(i)
                    newRimepop(i, j) = Best_rime(1, j);
                end

            end

        end

        for i = 1:SearchAgents_no
            %Boundary absorption
            Flag4ub = newRimepop(i, :) > ub;
            Flag4lb = newRimepop(i, :) < lb;
            newRimepop(i, :) = (newRimepop(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
            %newRime_rates(1,i)=fobj(newRimepop(i,:));
            if Time < Iteration
                [newRime_rates(1, i), FEs, newRimepop(i, :)] = sigle_evaluation(newRimepop(i, :), dim, thdim, fname, Pxy, FEs);

                %Positive greedy selection mechanism
            else
                break
            end

            if newRime_rates(1, i) > Rime_rates(1, i)
                Rime_rates(1, i) = newRime_rates(1, i);
                Rimepop(i, :) = newRimepop(i, :);

                %%
                if newRime_rates(1, i) > Best_rime_rate
                    Best_rime_rate = newRime_rates(1, i);
                    Best_rime = newRimepop(i, :);
                end

            end

        end

        Convergence_curve(Time) = Best_rime_rate;

        if Fbest < Best_rime_rate
            FE = FEs;
            iter = Time;
        end

        Fbest = Best_rime_rate;
        Lbest = Best_rime;

        Time = Time + 1;
    end
