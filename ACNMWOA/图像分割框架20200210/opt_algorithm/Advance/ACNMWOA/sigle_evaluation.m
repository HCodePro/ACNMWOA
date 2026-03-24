function [Fitness, fes, X] = sigle_evaluation(X, dim, thdim, fname, Pxy, fes)
    %SIGLE_CALCULATE 此处显示有关此函数的摘要
    %   此处显示详细说明
    %         Vtemp1=V(1,1:dim);
    %         Vtemp1(1,:)=sort(Vtemp1(1,:));
    %         Vtemp2=V(1,(dim+1):thdim);
    %         Vtemp2(1,:)=sort(Vtemp2(1,:));
    %         V=[Vtemp1 Vtemp2] ;
    %         FitnessV=feval(fname,V(1,:),Pxy);
    %         fes = fes + 1;
    %% 生成(1，dim)的随机矩阵
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
    fes = fes + 1;
end
