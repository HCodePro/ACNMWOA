function [best_pos,Convergence_curve]=SFOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
%%Hu, R., et al. (2017). "A short-term power load forecasting model based on the generalized regression neural network with decreasing step fruit fly optimization algorithm." Neurocomputing 221: 24-31.
tic
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs=0;
it=1;

%%Convergence_curve=zeros(1,Max_iter);
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


%% 果蝇寻优开始
%利用嗅觉寻找食物
for i = 1:sizepop
    tempx = ax*rand()-bx;
    tempy = ay*rand()-bx;
    for j = 1:dim
        X(i,j) = X_axis(1,j) + tempx;
        Y(i,j) = Y_axis(1,j) + tempy;
        D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
        S(i,j) = 1./D(i,j);
    end
    FEs = FEs + 1;
    Smell(i) = fobj(S(i,1:dim));
end

[bestSmell, bestindex] = min(Smell);
X_axis = X(bestindex,:);
Y_axis = Y(bestindex,:);
best_pos = S(bestindex,:);
bestCVaccuarcy = bestSmell;

%% 果蠅迭代尋優

 while FEs<MaxFEs  
    for i = 1:sizepop
        tempx = ax*rand()-bx;
        tempy = ay*rand()-by;
        for j = 1:dim
            X(i,j) = X_axis(1,j) + tempx*(1-(1/(1+exp(6-12*(FEs/MaxFEs)))));
            Y(i,j) = Y_axis(1,j) + tempy*(1-(1/(1+exp(6-12*(FEs/MaxFEs)))));
            D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
            S(i,j) = 1./D(i,j);
        end   
        FEs = FEs+1;
        Smell(i) = fobj(S(i,1:dim));
    end

	[bestSmell, bestindex] = min(Smell);
	if bestSmell < bestCVaccuarcy
        X_axis = X(bestindex,:);
        Y_axis = Y(bestindex,:);
        best_pos=S(bestindex,:);
        bestCVaccuarcy = bestSmell;
	end
     Convergence_curve(it)=(bestCVaccuarcy);
     it=it+1;
 end
toc
end

