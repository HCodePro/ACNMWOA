function [best_pos,Convergence_curve]=CIFOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
tic
%%An improved chaotic fruit fly optimization based on a mutation strategy for simultaneous feature selection and parameter optimization for SVM and its applications
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs=0;
it=1;
mr = 0.8;
u = (max(ub)-min(lb))/sizepop;
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
X_axis = sx*rand(1,dim);
Y_axis = sy*rands(1,dim);
X = zeros(sizepop, dim);
Y = zeros(sizepop, dim);
D = ones(sizepop, dim);
S = zeros(sizepop, dim);

C = zeros(sizepop, dim);
C0 = zeros(1, dim);

Smell = zeros(sizepop, 1); 

for i = 1:dim
    X_axis(1,i) = min(lb) + (max(ub) - min(lb)) * rand();
end
C0 = (X_axis-min(lb))/(max(ub)-min(lb));

for i = 1:sizepop
    for j = 1:dim
        if i == 1
            C(i,j) = 4.*C0(j).*(1-C0(j));
            C(i+1,j) = 4*C(i,j)*(1-C(i,j));
        else
            C(i+1,j) = 4*C(i,j)*(1-C(i,j));
        end
        
    end
end
% for i = 1:sizepop
%     for j = 1:dim
%         C(i,j) = lb+(ub-lb)*C(i,j);
%     end
% end
for i = 1:sizepop
    for j = 1:dim
        S(i,j) = min(X_axis)+(max(X_axis)-min(X_axis))*C(i,j);
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
best_pos = S(bestindex,:);
bestCVaccuarcy = bestSmell;

%% 果蠅迭代尋優

 while FEs<MaxFEs
     
    C0 = (X_axis-min(lb))/(max(ub)-min(lb));
    for i = 1:sizepop
       for j = 1:dim
          if i == 1
              C(i,j) = 4.*C0(j).*(1-C0(j));
              C(i+1,j) = 4*C(i,j)*(1-C(i,j));
          else
             C(i+1,j) = 4*C(i,j)*(1-C(i,j));
          end
       end
    end

    for i = 1:sizepop
        for j = 1:dim
            if rand()<=mr
                S(i,j) = min(X_axis(1,j))+(max(X_axis(1,j))-min(X_axis(1,j)))*C(i,j);
            else
                S(i,j) = X_axis(1,j)+sign(rand()-0.5)*rand()*u;
            end
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
        best_pos=S(bestindex,:);
        bestCVaccuarcy = bestSmell;
    else
        t= ceil(dim*rand());
        X_axis(1,t) = min(X_axis(1,t))+(max(X_axis(1,t))-min(X_axis(1,t)))*rand();
	end
     Convergence_curve(it)=(bestCVaccuarcy);
     it=it+1;
 end
toc
end

