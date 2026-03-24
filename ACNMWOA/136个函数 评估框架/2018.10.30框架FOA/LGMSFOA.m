function [best_pos,Convergence_curve]=LGMSFOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
%%Shan, D., et al. (2013). "LGMS-FOA: An improved fruit fly optimization algorithm for solving optimization problems." Mathematical Problems in Engineering 2013.
tic
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs=0;
it=1;

%%Convergence_curve=zeros(1,Max_iter);
%% 初始果蝇群体位置

X = zeros(sizepop, dim);
S = zeros(sizepop, dim);
Smell = zeros(sizepop, 1);
X_axis = 1*rands(1,dim);
n = 0.005;
oumiga0 = 1;
alpha = 0.95;
X_axis(1,:) = n * (min(lb) + (max(ub) - min(lb)) * rand());
%% 果蝇寻优开始
for i = 1:sizepop
    for j = 1:dim
%         X_axis(1,j) = n * ((rand()*(ub - lb))-lb);
        X(i,j) = X_axis(1,j) + oumiga0 * (min(lb) + (max(ub) - min(lb)) * rand());
        S(i,j) = X(i,j);
    end
    FEs = FEs+1;
    Smell(i) = fobj(S(i,1:dim));
end

[bestSmell, bestindex] = min(Smell);
X_axis = X(bestindex,:);
best_pos = S(bestindex,:);
bestCVaccuarcy = bestSmell;

%% 果蠅迭代尋優
 while FEs<MaxFEs
    oumiga = oumiga0 * (alpha^FEs);
    for i = 1:sizepop
        for j = 1:dim           
           X(i,j) = X_axis(1,j) + oumiga * (min(lb) + (max(ub) - min(lb)) * rand());
           S(i,j) = X(i,j);
        end
        FEs = FEs+1;
        Smell(i) = fobj(S(i,1:dim));
    end

    [bestSmell, bestindex] = min(Smell);
    if bestSmell < bestCVaccuarcy
        X_axis = X(bestindex,:);
        best_pos=S(bestindex,:);
        bestCVaccuarcy = bestSmell;
    end
     Convergence_curve(it)=(bestCVaccuarcy);
     it=it+1;
end
toc
end

