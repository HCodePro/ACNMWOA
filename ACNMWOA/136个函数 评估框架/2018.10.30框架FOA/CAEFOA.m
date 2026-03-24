function [best_pos,Convergence_curve]=CAEFOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
tic
%%Ship motion prediction using dynamic seasonal RvSVR with phase space reconstruction and the chaos adaptive efficient FOA
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs=0;
it =1;
lamda = 0.1;   %adjustment coefficient of step length for update
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

X_new = zeros(2*sizepop, dim);
Y_new = zeros(2*sizepop, dim);
D_new = ones(sizepop, dim);
S_new = zeros(sizepop, dim);
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
X_axis = X(bestindex,:);
Y_axis = Y(bestindex,:);
best_pos = S(bestindex,:);
bestCVaccuarcy = bestSmell;

%% 果蠅迭代尋優

 while FEs<MaxFEs
    
    u = MaxFEs/(lamda*(MaxFEs+FEs+1));
    for i = 1:sizepop
        tempx = ax*rand() - bx;
        tempy = ay*rand() - by;
        for j = 1:dim
            X(i,j) = X_axis(1,j) + u*min(abs(X_axis(1,j)-max(X_axis(1,j))),abs(X_axis(1,j)-min(X_axis(1,j))))+u*tempx;
            Y(i,j) = Y_axis(1,j) + u*min(abs(Y_axis(1,j)-max(Y_axis(1,j))),abs(Y_axis(1,j)-min(Y_axis(1,j))))+u*tempy;
        end
    end
    if rand()<=0.001
        P = zeros(2, sizepop);
        P(1,:) = 0.7;
        P(2,:) = 0.7;
        for i = 1:sizepop
            P(1,i+1) = mod((P(1,i)+P(2,i)),1);
            P(2,i+1) = mod((P(1,i)+2*P(2,i)),1);
        end
        for i = 1:sizepop
            X_new(i,:) = min(X_axis(1,:))+(max(X_axis(1,:))-min(X_axis(1,:))).*P(1,i);
            Y_new(i,:) = min(Y_axis(1,:))+(max(Y_axis(1,:))-min(Y_axis(1,:))).*P(2,i);
        end
        for i = sizepop+1:2*sizepop
            X_new(i,:) = X(i-sizepop,:);
            Y_new(i,:) = Y(i-sizepop,:);
        end
        for i = 1:(2*sizepop)
            for j = 1:dim
            D_new(i,j) = (X_new(i,j).^2 + Y_new(i,j).^2).^0.5;
            S_new(i,j) = 1./D_new(i,j);
            end
        end
        for i = 1:2*sizepop
            for j = 1:i
                if fobj(S_new(i,:))<fobj(S_new(j,:))
                    temp = X_new(i,:);
                    X_new(i,:) = X_new(j,:);
                    X_new(j,:) = temp;
                    temp = Y_new(i,:);
                    Y_new(i,:) = Y_new(j,:);
                    Y_new(j,:) = temp;
                end
            end
        end
        for i = 1:sizepop
            X(i,:) = X_new(i,:);
            Y(i,:) = Y_new(i,:);
        end
    end
     for i = 1:sizepop
         for j = 1:dim
            D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
            S(i,j) = 1./D(i,j);
         end
        Flag4ub=S(i,:)>ub;
        Flag4lb=S(i,:)<lb;
        S(i,:)=(S(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;   
        FEs = FEs+1;
        Smell(i) = fobj(S(i,1:dim));
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
     it=it+1;
 end
toc
end

