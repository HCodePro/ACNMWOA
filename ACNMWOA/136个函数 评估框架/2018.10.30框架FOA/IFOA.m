function [best_pos,Convergence_curve]=IFOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
tic
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs = 0;
it = 1;
%%Convergence_curve=zeros(1,MaxFEs);
%% 初始果蝇群体位置
foa_option.sx = 1;
foa_option.sy = 1;
foa_option.ax = 20;
foa_option.ay = 20;
foa_option.bx = 10;
foa_option.by = 10;
sx = foa_option.sx;
sy = foa_option.sy;
ax = foa_option.ax;
ay = foa_option.ay;
bx = foa_option.bx;
by = foa_option.by;
X_axis = sx*rands(1,dim);
Y_axis = sy*rands(1,dim);
X = zeros(sizepop, dim);
Y = zeros(sizepop, dim);
D = ones(sizepop, dim);
S = zeros(sizepop, dim);
Smell = zeros(sizepop, 1);

lamdaMax = 0.5 * (max(ub) - min(lb));
lamdaMin = 10^-5;

%% 果蝇寻优开始
for i = 1:sizepop
    for j = 1:dim
        S(i,j) = min(lb) + (max(ub) - min(lb)) * rand();
    end

    FEs = FEs+1;
    Smell(i) = fobj(S(i,1:dim));
end

%***根據初始味道濃度值尋找初始極值
[bestSmell, bestindex] = min(Smell);

%利用視覺尋找夥伴聚集味道濃度最高之處
%做法是保留最佳值初始位置及初始味道濃度
X_axis = X(bestindex,:);
Y_axis = Y(bestindex,:);

best_pos = S(bestindex,:);
bestCVaccuarcy = bestSmell;

%% 果蠅迭代尋優
while FEs<MaxFEs
    lamda = lamdaMax * exp(log(lamdaMin/lamdaMax) * (FEs/MaxFEs));
    
    %迭代求X，Y的混沌映射变量，即初始化混沌种群
    for i = 1:sizepop
        tempx = ax*rand() - bx;
        tempy = ay*rand() - by;
        for j = 1:dim
            X(i,j) = X_axis(1,j) + tempx;
            Y(i,j) = Y_axis(1,j) + tempy;
            D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
            S(i,j) = 1./D(i,j);
            d = ceil(rand()*dim); %取1到dim的随机整数
            if j == d
                S(i,j) = S(i,j) + sign(rand() - 0.5) * lamda;
            end
            if S(i,j) > ub
                S(i,j) = max(ub);
            end
            if S(i,j) < lb
                S(i,j) = min(lb);
            end
        end

        if FEs<MaxFEs
            FEs = FEs+1;
            Flag4ub=S(i,:)>ub;
            Flag4lb=S(i,:)<lb;
            S(i,:)=(S(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            Smell(i) = fobj(S(i,1:dim));
        end
    end
    
    %***根據味道濃度值尋找極值
    [bestSmell, bestindex] = min(Smell);
    
    %***迭代保留最佳值位置與味道濃度
    if bestSmell < bestCVaccuarcy
        X_axis = X(bestindex,:);
        Y_axis = Y(bestindex,:);
        best_pos=S(bestindex,:);
        bestCVaccuarcy = bestSmell;
    end
    Convergence_curve(it)=(bestCVaccuarcy);
    it = it+1;
end
toc
end

