function [best_pos,Convergence_curve]=AMFOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
tic
%%Melt index prediction by least squares support vector machines with an adaptive mutation fruit fly optimization algorithm
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs=0;
it=1;
%Convergence_curve=zeros(1,Max_iter);
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


X_axism = zeros(sizepop, dim);
Y_axism = zeros(sizepop, dim);
D1 = ones(sizepop, dim);
S1 = zeros(sizepop, dim);

%% 果蝇寻优开始
%利用嗅觉寻找食物
for i = 1:sizepop
    tempx = ax*rand() - bx;
    tempy = ay*rand() - by;
    for j = 1:dim
        X(i,j) = X_axis(1,j) + tempx;
        Y(i,j) = Y_axis(1,j) + tempy;
        D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
        S(i,j) = 1./D(i,j);
    end
	Flag4ub=S(i,:)>ub;
    Flag4lb=S(i,:)<lb;
    S(i,:)=(S(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
    FEs=FEs+1;
    Smell(i) = fobj(S(i,1:dim));
end

%***根據初始味道濃度值尋找初始極值
[bestSmell, bestindex] = min(Smell);

%利用視覺尋找夥伴聚集味道濃度最高之處
%做法是保留最佳值初始位置及初始味道濃度
bestCVaccuarcy = bestSmell;
X_axis = X(bestindex,:);
Y_axis = Y(bestindex,:);
best_pos=S(bestindex,:);
Smelltot = 0;
for k = 1:sizepop
   Smelltot = Smelltot+Smell(k);
end
Smell_avg = Smelltot/sizepop;
var = 0;
for k = 1:sizepop
    var = var+(Smell(k)-Smell_avg)*(Smell(k)-Smell_avg);
end
m = ceil(dim*rand());
Smell_1 = zeros(m, 1); 
if var>0 && bestCVaccuarcy>0
    for i = 1:m
        for j = 1:dim
        X_axism(i,j) = X_axis(1,j);
        Y_axism(i,j) = Y_axis(1,j);
        X_axism(i,j) = X_axis(1,j)+10*rand();
        Y_axism(i,j) = Y_axis(1,j)+10*rand();
        D1(i,j) = (X_axism(i,j).^2+Y_axism(i,j).^2).^0.5;
        S1(i,j) = 1./D1(i,j);
        end
        FEs=FEs+1;
        Smell_1(i) = fobj(S1(i,:));
    end
    [bestSmell_1, bestindex_1] = min(Smell_1);
    if bestSmell_1<bestCVaccuarcy
        bestCVaccuarcy = bestSmell_1;
        X_axis = X_axism(bestindex_1,:);
        Y_axis = Y_axism(bestindex_1,:);
        best_pos=S(bestindex,:);
    end
end
%% 果蠅迭代尋優

 while FEs<MaxFEs
    for i = 1:sizepop
        tempx = ax*rand() - bx;
        tempy = ay*rand() - by;
        for j = 1:dim
            X(i,j) = X_axis(1,j) + tempx;
            Y(i,j) = Y_axis(1,j) + tempy;
            D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
            S(i,j) = 1./D(i,j);
        end
        if FEs<MaxFEs
            Flag4ub=S(i,:)>ub;
            Flag4lb=S(i,:)<lb;
            S(i,:)=(S(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
            FEs=FEs+1;
            Smell(i) = fobj(S(i,1:dim));
        end
    end
    
	%***根據味道濃度值尋找極值
	[bestSmell, bestindex] = min(Smell);
  
	%***迭代保留最佳值位置與味道濃度
    if bestSmell < bestCVaccuarcy
        X_axis = X(bestindex,:);
        Y_axis = Y(bestindex,:);
        bestCVaccuarcy = bestSmell;
    end
    
    

Smelltot = 0;
for k = 1:sizepop
   Smelltot = Smelltot+Smell(k);
end
Smell_avg = Smelltot/sizepop;
var = 0;
for k = 1:sizepop
    var = var+(Smell(k)-Smell_avg)*(Smell(k)-Smell_avg);
end
m = ceil(dim*rand());
Smell_1 = zeros(m, 1); 
if var>0 && bestCVaccuarcy>0
    for i = 1:m
        for j = 1:dim
        X_axism(i,j) = X_axis(1,j);
        Y_axism(i,j) = Y_axis(1,j);
        X_axism(i,j) = X_axis(1,j)+10*rand();
        Y_axism(i,j) = Y_axis(1,j)+10*rand();
        D1(i,j) = (X_axism(i,j).^2+Y_axism(i,j).^2).^0.5;
        S1(i,j) = 1./D1(i,j);
        end
        FEs = FEs+1;
        Smell_1(i) = fobj(S1(i,:));
    end
    [bestSmell_1, bestindex_1] = min(Smell_1);
    if bestSmell_1<bestCVaccuarcy
        bestCVaccuarcy = bestSmell_1;
        X_axis = X_axism(bestindex_1,:);
        Y_axis = Y_axism(bestindex_1,:);
    end
end
     Convergence_curve(it)=(bestCVaccuarcy);
     it=it+1;
 end
toc
end

