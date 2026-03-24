% function [BestSol,Convergence_curve]=DE( SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
function [Fbest,Lbest,FE,MaxFEs,Convergence_curve,iter]=IterationEDE_100(fname,SearchAgents_no,thdim,lb,ub,MaxFEs,Pxy,Iteration)
% 参数向量 parameters [n,N_iteration,beta_min,beta_max,pCR]
% n为种群规模，N_iteration为迭代次数
% beta_min 缩放因子下界 Lower Bound of Scaling Factor
% beta_max=0.8; % 缩放因子上界 Upper Bound of Scaling Factor
% pCR 交叉概率 Crossover Probability
% 要求输入数据为列向量（矩阵）
%%
dim=thdim/2;
Fbest=-inf;

k=SearchAgents_no/2;
trial=zeros(1,SearchAgents_no/2);
limit=100;
%%
para=[SearchAgents_no,MaxFEs,0.2,0.8,0.2];
%% 差分进化（DE）算法
nPop=para(1); % 种群规模 Population Size
MaxFEs=para(2); % 最大迭代次数Maximum Number of Iterations
nVar=dim; % 自变量维数，此例需要优化两个参数c和g Number of Decision Variables
VarSize=[1,2*dim]; % 决策变量矩阵大小 Decision Variables Matrix Size
beta_min=para(3); % 缩放因子下界 Lower Bound of Scaling Factor
beta_max=para(4); % 缩放因子上界 Upper Bound of Scaling Factor
pCR=para(5); %  交叉概率 Crossover Probability
lb=ones(1,2*dim).*lb; % 参数取值下界
ub=ones(1,2*dim).*ub; % 参数取值上界
%% 初始化 Initialization
FEs=0;
N=para(1)/2;
bestPositions=zeros(1,2*dim);
Destination_fitness=-inf;
AllFitness = -inf*ones(N,1);
X=random_initialization(nPop,dim,ub,lb);
for i=1:nPop/2 
    [AllFitness(i),FEs,X(i,:)] =sigle_evaluation(X(i,:),dim,thdim,fname,Pxy,FEs);      
%     AllFitness(i)=fobj(X(i,:));
%     FEs=FEs+1;
    if AllFitness(i) > Destination_fitness 
        bestPositions=X(i,:); 
        Destination_fitness=AllFitness(i);
    end    
end
BestCost=zeros(MaxFEs,1); 
Convergence_curve=[];
it=1;
%% 主循环 DE Main Loop
while  FEs < MaxFEs || it<=Iteration
    for i=1:nPop/2 % 遍历每个个体
        x=X(i,:); % 提取个体位置
        % 随机选择三个个体以备变异使用
        A=randperm(nPop/2); % 个体顺序重新随机排列
        A(A==i)=[]; % 当前个体所排位置腾空（产生变异中间体时当前个体不参与）
        a=A(1);
        b=A(2);
        c=A(3);
        % 变异操作 Mutation
        beta=unifrnd(beta_min,beta_max,VarSize); % 随机产生缩放因子
        y=X(a,:)+beta.*(X(b,:)-X(c,:)); % 产生中间体
        % 防止中间体越界
        y=max(y,lb);
		y=min(y,ub);
        % 交叉操作 Crossover
        z=zeros(size(x)); % 初始化一个新个体
        j0=randi([1,numel(x)]); % 产生一个伪随机数，即选取待交换维度编号？？？
        for j=1:numel(x) % 遍历每个维度
            if j==j0 || rand<=pCR % 如果当前维度是待交换维度或者随机概率小于交叉概率
                z(j)=y(j); % 新个体当前维度值等于中间体对应维度值
            else
                z(j)=x(j); % 新个体当前维度值等于当前个体对应维度值
            end
        end
        
         [Fz,FEs,z] =sigle_evaluation(z,dim,thdim,fname,Pxy,FEs);   
%         Fz=fobj(z); % 新个体目标函数值
%         FEs=FEs+1;
        if Fz>AllFitness(i) % 如果新个体优于当前个体
            X(i,:)=z; % 更新当前个体
            AllFitness(i)=Fz;
            trial(i)=0;
            if Fz > Destination_fitness % 如果当前个体（更新后的）优于最优个体
               bestPositions = z; % 更新最优个体
               Destination_fitness = Fz;
            end
        else
            trial(i)=trial(i)+1;
        end
    end  
    %% 嵌入ABC更新框架，但不使用其更新方程 
    i=1;
    t=0;
    while t<nPop/2         
        t=t+1;          
        Coef=FEs/MaxFEs;          %% 可改进
        r1=rand();                         
        Beta=exp(r1*((MaxFEs-FEs+1)/MaxFEs))*(sin(2*pi*r1));     %% 可改进
        sol=ones(1,2*dim);
        if t==1
            Temp=bestPositions;
        else
            Temp=X(t-1,:);
        end
        if  Coef > rand %shrink surround strategy                                                     
            sol=bestPositions+rand(1,2*dim).*(Temp-X(t,:))+Beta*(bestPositions-X(t,:)); 
        else
            IndivRand=rand(1,2*dim).*(ub-lb)+lb;                                
            sol=IndivRand+rand(1,2*dim).*(Temp-X(t,:))+Beta*(IndivRand-X(t,:));  
        end  
        
        ind = find(sol<lb);
        sol(ind) = lb(ind);
        ind = find(sol>ub);
        sol(ind) = ub(ind);
        [F_sol,FEs,sol] =sigle_evaluation(sol,dim,thdim,fname,Pxy,FEs);           
%         F_sol=fobj(sol);
%         FEs = FEs +1;
        if F_sol > AllFitness(i)
            X(i,:)=sol;
            AllFitness(i)=F_sol;
            if F_sol > Destination_fitness
                bestPositions = sol;
                Destination_fitness = F_sol;
            end
            trial(i)=0;
        else
            trial(i)=trial(i)+1; 
        end
        i=i+1;
        if (i==(nPop/2)+1)
            i=1;     
        end
    end
    %% 检验是否还具有价值
    ind=find(trial==max(trial));
    ind=ind(end);
    if  trial(ind) > limit
        trial(ind)=0;
        sol=(ub-lb).*rand(1,2*dim)+lb;
        [F_sol,FEs,sol] =sigle_evaluation(sol,dim,thdim,fname,Pxy,FEs);             
%         F_sol=fobj(sol);
%         FEs = FEs + 1;
        X(ind,:)=sol;
        AllFitness(ind)=F_sol;
    end
    %% 将适应度值排序，选取前k个解
    [X_fitnessed, SortOrder]=sort(AllFitness,'descend');
    X=X(SortOrder,:);
    AllFitness=X_fitnessed(1:N);
    X_Elite_Temp=X(1:k,:);     
    %%  Elite Levy Spreading Strategy
    eta=0.4*exp(-(2*FEs/MaxFEs)^2);
    newX=ones(1,2*dim);
    for i=1:k
        kk=rand(1,2*dim)>eta;
        Rand_k=floor(rand*k+1);
        Rand_N=floor(rand*N+1);  %% 促进交流
        Float=(2*rand()-1)*(1-(FEs/MaxFEs));
        for j=1:2*dim
            step=Float*(X_Elite_Temp(Rand_k,j)-X(Rand_N,j));
            newX(j)=X(i,j)+step*kk(j)*Levy(1);
            if newX(j)>ub(j) || newX(j)<lb(j)
                newX(j) = eta*bestPositions(j)+(1-eta)*(rand()*(ub(j)-lb(j)));
            end
        end
        
        [F_newX,FEs,newX] =sigle_evaluation(newX,dim,thdim,fname,Pxy,FEs);           
%         F_newX = fobj(newX);
%         FEs=FEs+1;
        if F_newX > AllFitness(i)
            X(i,:)=newX;
            AllFitness(i) = F_newX;
            if F_newX > Destination_fitness
                bestPositions = newX;
                Destination_fitness = F_newX;
            end
        end
    end
    k=1+floor((1-FEs/MaxFEs)*(nPop/2-1));
    %%
    BestCost(it)=Destination_fitness;  
    Convergence_curve(it)=Destination_fitness;
    
    if Fbest<Destination_fitness
    FE=FEs;
    iter=it;
    end
    Fbest=Destination_fitness;
    Lbest=bestPositions;
    it=it+1;
end
Positionbest=bestPositions;
bestCVaccuarcy=Destination_fitness;
end

function o=Levy(d)
beta=1.5;
%Eq. (3.10)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,2*d)*sigma;
v=randn(1,2*d);
step=u./abs(v).^(1/beta);
% Eq. (3.9)
o=step;
end