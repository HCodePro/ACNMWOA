function [best_pos,Convergence_curve]=FOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
tic
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs=0;
it=1;
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
    for i = 1:sizepop
        tempx = ax*rand() - bx;
        tempy = ay*rand() - by;
        for j = 1:dim
            X(i,j) = X_axis(1,j) + tempx;
            Y(i,j) = Y_axis(1,j) + tempy;
            D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
            S(i,j) = 1./D(i,j);
        end
        
%         Moth_pos_m_gaus=S(i,:)*(1+Levy1(1));
%         FEs=FEs+2;
%         
%         Moth_fitness_m_gaus=fobj(Moth_pos_m_gaus);
%         Moth_fitness_s=fobj(S(i,:));
%         %         Moth_fitness_f=fobj(Moth_pos(i,:));
%         
%         Moth_fitness_comb=[Moth_fitness_m_gaus,Moth_fitness_s];
%         [~,m]=min(Moth_fitness_comb);
%         if m==1
%             S(i,:)=Moth_pos_m_gaus;
%         else
%             S(i,:)=S(i,:);
%         end

        
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
        best_pos=S(bestindex,:);
        bestCVaccuarcy = bestSmell;
	end
     Convergence_curve(it)=(bestCVaccuarcy);
     it=it+1;
 end
toc
end


% function o=Levy1(d)
% % d=1;
% beta=3/2;
% %Eq. (3.10)
% sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
% u=randn(1,d)*sigma;
% v=randn(1,d);
% step=u./abs(v).^(1/beta);
% 
% % Eq. (3.9)
% o=step;
% 
% end

