% function [BestSol,Convergence_curve]=DE( SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
function [Fbest,Lbest,FE,MaxFEs,Convergence_curve,iter]=EvolutionGLSMA(fname,N,thdim,lb,ub,MaxFEs,Pxy,Iteration)
%% Initialize position
dim=thdim/2;
Fbest=-inf;

bestPositions=zeros(1,2*dim);
Destination_fitness=-inf;
AllFitness = -inf*ones(N,1);
weight = ones(N,2*dim);
%% Initialize the set of random solutions
Convergence_curve=[];
it=1;
FEs=1;
X=random_initialization(N,dim,ub,lb);
X1 = X;
z=0.03; % 参数
for i=1:N
    Flag4ub=X(i,:)>ub;
    Flag4lb=X(i,:)<lb;
    X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    FEs = FEs+1;
%     AllFitness(i) = fobj(X(i,:));
    [AllFitness(i),FEs,X(i,:)] =sigle_evaluation(X(i,:),dim,thdim,fname,Pxy,FEs);
end
[SmellOrder,SmellIndex] = sort(AllFitness,'descend');
worstFitness = SmellOrder(size(X,1));%初设最后一个是最差的
bestFitness = SmellOrder(1);%初设第一个是最好的
S=bestFitness-worstFitness+eps;
if bestFitness > Destination_fitness
        bestPositions=X1(SmellIndex(1),:);
        Destination_fitness = bestFitness;
end

count = 0;%计数
%% Main loop
while  FEs < MaxFEs || it<=Iteration
    %% SMA更新操作
    for i=1:N
        for j=1:2*dim       
            if i<=(N/2)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    flag = true;
    %更新最佳适应度值和最佳位置
    if bestFitness > Destination_fitness
        bestPositions=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
        flag = false;
    end
    
    %%--如果一直没有更新，就使用莱维飞行机制
    if flag 
        count = count + 1;
    else
        count = 0;
    end
    
    if count == 10
        %执行某个搜索
        X_levy = X(i,:)+rand()*sign(rand()-1/2).*Levy1(2*dim);%%莱维飞行
        Flag4ub=X_levy>ub;
        Flag4lb=X_levy<lb;
        X(i,:)=X_levy.*(~(Flag4ub+Flag4lb))+ub.*Flag4ub+lb.*Flag4lb;
    end
    
    a = atanh(-(FEs/MaxFEs)+1);  
    b = 1-FEs/MaxFEs;
    XX = zeros(1,2*dim);
    %%
    for i=1:N
        if rand<z     %Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness)); %Eq.(2.2)
            vb = unifrnd(-a,a,1,2*dim);  %Eq.(2.3)
            vc = unifrnd(-b,b,1,2*dim);
            for j=1:2*dim
                r = rand();
                A = randi([1,N]);    % 从种群中随机挑选两个个体
                B = randi([1,N]);
                if r<p          %Eq.(2.1)
                    X(i,j) = bestPositions(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                else
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
         %%%%%%%%%%%%%%%%%%%%%%%高斯变异机制 
         X_gaus=X(i,:).*(1+randn(1));    %randn()是一种产生标准正态分布的随机数或矩阵的函数，                                                  
         Flag4ub=X_gaus>ub;
         Flag4lb=X_gaus<lb;
         X_gaus=(X_gaus.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;         
         [X_fitness_gaus,FEs,X_gaus] =sigle_evaluation(X_gaus,dim,thdim,fname,Pxy,FEs);

         Flag4ub=X(i,:)>ub;
         Flag4lb=X(i,:)<lb;
         X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;    
         [X_fitness_s,FEs,X(i,:)] =sigle_evaluation(X(i,:),dim,thdim,fname,Pxy,FEs);         
         X_fitness_comb=[X_fitness_gaus,X_fitness_s];
         [~,m]=max(X_fitness_comb);
         if m==1
             X(i,:)=X_gaus;
         else
             X(i,:)=X(i,:);
         end  
    end
    
    for i=1:N
        [temp_fit,FEs,X(i,:)] =sigle_evaluation(X(i,:),dim,thdim,fname,Pxy,FEs);
        if AllFitness(i) < temp_fit
			AllFitness(i) = temp_fit;
			X1(i,:) = X(i,:);
        end
    end   
	[SmellOrder,SmellIndex] = sort(AllFitness,'descend');
    worstFitness = SmellOrder(size(X,1));
    bestFitness = SmellOrder(1);
    S=bestFitness-worstFitness+eps;
    
    if bestFitness > Destination_fitness
        bestPositions=X1(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end    
    %% 
    Convergence_curve(it)=Destination_fitness;
    if Fbest<Destination_fitness
        FE=FEs;
        iter=it;
    end
    Fbest=Destination_fitness;
    Lbest=bestPositions;
    it=it+1;
end
end
%莱维飞行
function o=Levy1(d)
    beta=1.5;  %beta设置为1.5，与LWOA23相同
    %Eq. (3.10)
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;
    v=randn(1,d);
    step=u./abs(v).^(1/beta);

    % Eq. (3.9)
    o=step;
end