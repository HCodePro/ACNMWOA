% function [Destination_position,Convergence_curve]=m_SCA(N,MaxFEs,lb,ub,2*dim,fobj)
function [Fbest, Lbest, FE, MaxFEs, Convergence_curve, iter] = m_SCA(fname, N, thdim, lb, ub, MaxFEs, Pxy, Iteration)

    if Iteration ~= 0
        MaxFEs = Iteration;
    end

    dim = thdim / 2;
    Fbest = -inf;

    FEs = 0;
    %Initialize the set of random solutions
    % X=initialization(N,2*dim,ub,lb);
    X = random_initialization(N, dim, ub, lb);
    Destination_position = zeros(1, 2 * dim);
    Destination_fitness = -inf;

    Convergence_curve = [];
    Objective_values = zeros(1, size(X, 1));

    % Calculate the fitness of the first set and find the best one
    for i = 1:size(X, 1)
        %     Objective_values(1,i)=fobj(X(i,:));
        %     FEs=FEs+1;

        [Objective_values(1, i), FEs, X(i, :)] = sigle_evaluation_m_SCA(X(i, :), dim, thdim, fname, Pxy, FEs, Iteration);

        if i == 1
            Destination_position = X(i, :);
            Destination_fitness = Objective_values(1, i);
        elseif Objective_values(1, i) > Destination_fitness
            Destination_position = X(i, :);
            Destination_fitness = Objective_values(1, i);
        end

    end

    %Initialize the algorithm parameters – N (size of population of solutions), T (maximum number of iterations as termination criteria), A , perturbation or jumping rateJ R and self-adaptation rate S R
    JR = 0.1;
    %Main loop
    t = 1; % start from the second iteration since the first iteration was dedicated to calculating the fitness

    while FEs <= MaxFEs

        if Iteration ~= 0
            FEs = t;
        end

        %     display('m_SCA');
        %     display(t);

        % Eq. (3.4)
        a = 2;
        r1 = a - a * (FEs / MaxFEs); % r1 decreases linearly from a to 0

        if (rand < JR)
            [P, FEs] = OBL(X, N, 2 * dim, ub, lb, 0, FEs, thdim, fname, Pxy, Iteration);

            X = P;

            for i = 1:size(X, 1)
                % Check if solutions go outside the search spaceand bring them back
                Flag4ub = X(i, :) > ub;
                Flag4lb = X(i, :) < lb;
                X(i, :) = (X(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

                % Calculate the objective values
                %             Objective_values(1,i)=fobj(X(i,:));
                %FEs=FEs+1;
                [Objective_values(1, i), FEs, X(i, :)] = sigle_evaluation_m_SCA(X(i, :), dim, thdim, fname, Pxy, FEs, Iteration);

                % Update the destination if there is a better solution
                if Objective_values(1, i) > Destination_fitness
                    Destination_position = X(i, :);
                    Destination_fitness = Objective_values(1, i);
                end

            end

        else
            % Update the position of solutions with respect to destination
            for i = 1:size(X, 1) % in i-th solution

                for j = 1:size(X, 2) % in j-th 2*dimension

                    % Update r2, r3, and r4 for Eq. (3.3)
                    r2 = (2 * pi) * rand();
                    r3 = 2 * rand;
                    r4 = rand();
                    SR = rand();

                    % Eq. (3.3)
                    if r4 < 0.5
                        % Eq. (3.1)
                        X(i, j) = X(i, j) + (r1 * sin(r2) * abs(r3 * Destination_position(j) - X(i, j))) + SR * (Destination_position(j) - X(i, j));
                    else
                        % Eq. (3.2)
                        X(i, j) = X(i, j) + (r1 * cos(r2) * abs(r3 * Destination_position(j) - X(i, j))) + SR * (Destination_position(j) - X(i, j));
                    end

                end

            end

            for i = 1:size(X, 1)
                % Check if solutions go outside the search spaceand bring them back
                Flag4ub = X(i, :) > ub;
                Flag4lb = X(i, :) < lb;
                X(i, :) = (X(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

                % Calculate the objective values
                %             Objective_values(1,i)=fobj(X(i,:));
                %             FEs=FEs+1;
                [Objective_values(1, i), FEs, X(i, :)] = sigle_evaluation_m_SCA(X(i, :), dim, thdim, fname, Pxy, FEs, Iteration);

                % Update the destination if there is a better solution
                if Objective_values(1, i) > Destination_fitness
                    Destination_position = X(i, :);
                    Destination_fitness = Objective_values(1, i);
                end

            end

        end

        Convergence_curve(t) = Destination_fitness;
        %    display(['At iteration ', num2str(t), ' the best universes fitness is ', num2str(Convergence_curve(t))]);
        if Fbest < Destination_fitness
            FE = FEs;
            iter = t;
        end

        Fbest = Destination_fitness;
        Lbest = Destination_position;
        % Increase the iteration counter
        t = t + 1;
    end

end

function [Fitness, fes, X] = sigle_evaluation_m_SCA(X, dim, thdim, fname, Pxy, fes, Iteration)
    %SIGLE_CALCULATE 此处显示有关此函数的摘要
    %   此处显示详细说明
    %         Vtemp1=V(1,1:dim);
    %         Vtemp1(1,:)=sort(Vtemp1(1,:));
    %         Vtemp2=V(1,(dim+1):thdim);
    %         Vtemp2(1,:)=sort(Vtemp2(1,:));
    %         V=[Vtemp1 Vtemp2] ;
    %         FitnessV=feval(fname,V(1,:),Pxy);
    %         fes = fes + 1;
    X = round(X);
    Xtemp1 = X(:, 1:dim);

    for si = 1:size(Xtemp1, 1)
        Xtemp1(si, :) = sort(Xtemp1(si, :));
    end

    Xtemp2 = X(:, (dim + 1):thdim);

    for si = 1:size(Xtemp2, 1)
        Xtemp2(si, :) = sort(Xtemp2(si, :));
    end

    X = [Xtemp1 Xtemp2];
    Fitness = feval(fname, X(1, :), Pxy);

    if Iteration == 0
        fes = fes + 1;
    end

end
