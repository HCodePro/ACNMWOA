function [P, FEs] = EOBL(X, dim, ub, lb, thdim, fname, Pxy, FEs)
    %EOBL 此处显示有关此函数的摘要
    %   此处显示详细说明

    ub_oppo = max(X);
    lb_oppo = min(X);

    for i = 1:size(X, 1)
        [fitness(1, i), FEs, X(i, :)] = sigle_evaluation(X(i, :), dim, thdim, fname, Pxy, FEs);
    end

    for i = 1:size(X, 1)
        GOX(i, :) = rand() * (ub_oppo + lb_oppo) - X(i, :);

        for j = 1:2 * dim

            if GOX(i, j) > ub_oppo(j) || GOX(i, j) < lb_oppo(j)
                GOX(i, j) = rand() * (ub_oppo(j) - lb_oppo(j)) + lb_oppo(j);
            end

        end

        [Gfitness(1, i), FEs, GOX(i, :)] = sigle_evaluation(GOX(i, :), dim, thdim, fname, Pxy, FEs);
    end

    allX = [X; GOX];
    allfitness = [fitness, Gfitness];
    [~, index] = sort(allfitness, 'descend');
    P = allX(index, :);
    P = P(1:size(X, 1), :);

end
