function [best_pos,Convergence_curve]=FOASR_GS_Chaotic(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
tic
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs=0;
it =1;
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

S_temp = zeros(sizepop, dim);
SCauss = zeros(sizepop, dim);
ChaosVec = zeros(MaxFEs,1);
B = zeros(MaxFEs,1);
ChaosVec(1)=0.95;


% Iterative map
a=0.7;
for i=1:MaxFEs
    ChaosVec(i+1)=sin((a*pi)/ChaosVec(i));
    B(i)=((ChaosVec(i)+1)*1)/2;
    %     G(i)=lb+(ub-lb)*B(i);
end
ChaosVec = B;

% %Sine map
% for i=1:Max_iter
%     ChaosVec(i+1) = sin(pi*ChaosVec(i));
%     %      G(i)=(x(i))*Value;
% end

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
    S_temp(i,:) = S(i,:);
	Flag4ub=S(i,:)>ub;
    Flag4lb=S(i,:)<lb;
    S(i,:)=(S(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;   
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
    for i = 1:sizepop
        tempx = ax*rand() - bx;
        tempy = ay*rand() - by;
        change_step = 0.8;
        while change_step<=1
            change_step = change_step+0.1;
            for j = 1:dim
              X(i,j) = X_axis(1,j) + tempx*change_step;
              Y(i,j) = Y_axis(1,j) + tempy*change_step;
              D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
              S_temp(i,j) = 1./D(i,j);
            end
            if S_temp(i,:) > ub
            S_temp(i,:) = ub;
            end
            if S_temp(i,:) < lb
            S_temp(i,:) = lb;
            end
            if FEs<MaxFEs
                FEs = FEs+2;
                fitness_st = fobj(S(i,:));
                fitness_temp = fobj(S_temp(i,:));
                fitness_comp = [fitness_temp,fitness_st];
                [~,p] = min(fitness_comp);
                if p==1
                S(i,:) = S_temp(i,:);
                else
                S(i,:) = S(i,:);
                end 
            end
        end
        
        if FEs<MaxFEs
            FEs = FEs+4;
            SChaos = S(i,:) + ChaosVec(FEs)* (S(i,:) - best_pos(1,:));
            SChaos_fitness = fobj(SChaos);
            Moth_pos_m_gaus=S(i,:)*(1+randn(1));
            Moth_fitness_m_gaus=fobj(Moth_pos_m_gaus);
            Moth_fitness_s=fobj(S(i,:));
            Moth_fitness_comb=[Moth_fitness_m_gaus,Moth_fitness_s];
            [~,m]=min(Moth_fitness_comb);
            if m==1    
                cfoa_fitness_comb=[SChaos_fitness,Moth_fitness_m_gaus];
                [~,w]=min(cfoa_fitness_comb);
                if w==1
                    S(i,:)=SChaos;
                else
                    S(i,:)=Moth_pos_m_gaus;
                end
            else
            S(i,:)=S(i,:); 
            end     
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
     it=it+1;
 end
toc
end

