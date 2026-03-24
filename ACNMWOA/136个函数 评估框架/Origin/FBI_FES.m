% Composed and written by Nguyen Ngoc Mai - PiM Lab - National Taiwan University of Science and Technology
% Supervisor: Distinguished Professor Jui-Sheng Chou
% Citation:
% Jui-Sheng Chou, Ngoc-Mai Nguyen, FBI Inspired Meta-Optimization, Applied Soft Computing Journal (2020) 106339, doi: https://doi.org/10.1016/j.asoc.2020.106339.
% ---------------------------------------------------------------------------------------------------------
function [X, Fit_store] = FBI_FES(NP, GEN, LB, UB, D, fobj)
    %% Initialized population
    ObjVal = zeros(1, NP);
    LB = ones(1, D) .* LB;
    UB = ones(1, D) .* UB;
    FES = 0;
    Pop = initialization(NP, D, UB, LB);

    for i = 1:NP
        ObjVal(1, i) = fobj(Pop(i, :));
        FES = FES + 1;
    end

    %% Memorize the best solution
    iBest = find(ObjVal == min(ObjVal)); %find：查找非零元素的索引  这里是查找特定值的索引
    iBest = iBest(end);
    GlobalMin = ObjVal(iBest);
    Xbest = Pop(iBest, :);
    g = 1;
    Fit_store = [];
    %% Optimization Cycle
    PopA = Pop; ObjValA = ObjVal;
    PopB = Pop; ObjValB = ObjVal;

    while FES <= GEN
        %% Investigation team - team A
        %% Step A1
        for i = 1:NP
            Change = fix(rand * D) + 1; %fix:选取向0靠近的整数
            nb1 = floor(rand * NP) + 1;

            while (nb1 == i)
                nb1 = floor(rand * NP) + 1;
            end

            nb2 = floor(rand * NP) + 1;

            while (nb1 == nb2 || nb2 == i)
                nb2 = floor(rand * NP) + 1;
            end

            solA = PopA(i, :);
            solA(Change) = PopA(i, Change) + (PopA(i, Change) - (PopA(nb1, Change) + PopA(nb2, Change)) / 2) * (rand - 0.5) * 2; %Eq.(2) in FBI Inspired Meta-Optimization

            for Change = 1:D

                if (solA(Change) <= LB(Change)) || (solA(Change) >= UB(Change))
                    solA(Change) = LB(Change) + (UB(Change) - LB(Change)) * rand();
                end

            end

            f_a = fobj(solA);
            FES = FES + 1;

            if f_a <= ObjValA(i)
                PopA(i, :) = solA;
                ObjValA(i) = f_a;

                if f_a <= GlobalMin
                    Xbest = solA;
                    GlobalMin = f_a;
                end

            end

        end

        %% Step A2
        if min(ObjValA) < max(ObjValA)
            prob = probability(ObjValA);

            for i = 1:NP

                if (rand > prob(i))
                    r(1) = floor(rand() * NP) + 1;

                    while r(1) == i
                        r(1) = floor(rand() * NP) + 1;
                    end

                    r(2) = floor(rand() * NP) + 1;

                    while (r(2) == r(1)) || (r(2) == i)
                        r(2) = floor(rand() * NP) + 1;
                    end

                    r(3) = floor(rand() * NP) + 1;

                    while (r(3) == r(2)) || (r(3) == r(1)) || (r(3) == i)
                        r(3) = floor(rand() * NP) + 1;
                    end

                    solA = PopA(i, 1:D);
                    Rnd = floor(rand() * D) + 1;

                    for j = 1:D

                        if (rand() < rand()) || (Rnd == j)
                            solA(j) = Xbest(j) + PopA(r(1), j) + rand() * (PopA(r(2), j) - PopA(r(3), j)); %Eq.(5) in FBI Inspired Meta-Optimization
                        else
                            solA(j) = PopA(i, j);
                        end

                    end

                    for j = 1:D

                        if (solA(j) <= LB(j)) || (solA(j) >= UB(j))
                            solA(j) = LB(j) + (UB(j) - LB(j)) * rand();
                        end

                    end

                    f_a = fobj(solA);
                    FES = FES + 1;

                    if f_a <= ObjValA(i)
                        PopA(i, :) = solA;
                        ObjValA(i) = f_a;

                        if f_a <= GlobalMin
                            Xbest = solA;
                        end

                    end

                end

            end

        end

        %% Persuing team - team B
        %% Step B1
        for i = 1:NP
            SolB = PopB(i, 1:D);

            for j = 1:D
                SolB(j) = rand() * PopB(i, j) + rand() * (Xbest(j) - PopB(i, j)); %Eq.(6) in FBI Inspired Meta-Optimization

                if (SolB(j) < LB(j)) || (SolB(j) > UB(j))
                    SolB(j) = LB(j) + (UB(j) - LB(j)) * rand();
                end

            end

            f_b = fobj(SolB);
            FES = FES + 1;

            if f_b <= ObjValB(i)
                PopB(i, 1:D) = SolB;
                ObjValB(i) = f_b;

                if f_b <= GlobalMin
                    Xbest = SolB;
                    GlobalMin = f_b;
                end

            end

        end

        %% Step B2
        for i = 1:NP
            rr = floor(rand() * NP) + 1;

            while rr == i
                rr = floor(rand() * NP) + 1;
            end

            if ObjValB(i) > ObjValB(rr)
                SolB = PopB(rr, :) + rand(1, D) .* (PopB(rr, :) - PopB(i, :)) + rand(1, D) .* (Xbest - PopB(rr, :)); %Eq.(7) in FBI Inspired Meta-Optimization

                for j = 1:D

                    if (SolB(j) < LB(j)) || (SolB(j) > UB(j))
                        SolB(j) = LB(j) + (UB(j) - LB(j)) * rand();
                    end

                end

            else
                SolB = PopB(i, :) + rand(1, D) .* (PopB(i, :) - PopB(rr, :)) + rand(1, D) .* (Xbest - PopB(i, :)); %Eq.(8) in FBI Inspired Meta-Optimization

                for j = 1:D

                    if (SolB(j) < LB(j)) || (SolB(j) > UB(j))
                        SolB(j) = LB(j) + (UB(j) - LB(j)) * rand();
                    end

                end

            end

            f_b = fobj(SolB);
            FES = FES + 1;

            if f_b <= ObjValB(i)
                PopB(i, :) = SolB;
                ObjValB(i) = f_b;

                if f_b <= GlobalMin
                    Xbest = SolB;
                    GlobalMin = f_b;
                end

            end

        end

        %     if rem(g, ishow) == 0
        %         fprintf('Generation: %d. Best f: %f.\n', g, GlobalMin);
        %     end
        Fit_store(g) = GlobalMin;
        g = g + 1;
    end

    %% Result
    f = GlobalMin;
    X = Xbest;
end

% %% functions
% function f = rescale_matrix(X, LB, UB)
% [NP,D] = size(X);
% f = zeros(NP,D);
% for i = 1:D
%     f(:, i) = LB(i)*ones(NP,1) + (UB(i) - LB(i))*X(:,i);
% end
% end
%% functions
function prob = probability(fObjV)
    prob = (max(fObjV) - fObjV) / ((max(fObjV) - min(fObjV))); %Eq.(3) in FBI Inspired Meta-Optimization
end
