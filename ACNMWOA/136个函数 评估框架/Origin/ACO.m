% function [Leader_pos,Convergence_curve]=CC_ACO_R(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
function [Leader_pos,Convergence_curve]=ACO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)

%% Problem Definition
solution=[1 dim];             %初始化一个解决方法的大小


%% 初始化相关参数 
FEs=0;
it=1;
k=9;            % 档案袋大小k
m=SearchAgents_no;             % 新生成个数
q=0.6;
ibslo=1;     
p2=0.1;
%% 初始化
% 初始化一个个体
empty_individual.Position=[];
empty_individual.fitness=[];


% fitness_history=zeros(k,Max_iteration);
% position_history=zeros(k,Max_iteration,dim);
% Trajectories=zeros(k,Max_iteration);
% 初始化种群
pop=repmat(empty_individual,k,1);
for i=1:k
    pop(i).Position=unifrnd(lb,ub,solution);
    Flag4ub=pop(i).Position>ub;        %返回大于ub的逻辑值
    Flag4lb=pop(i).Position<lb;
    pop(i).Position=(pop(i).Position.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; 
    pop(i).fitness=fobj(pop(i).Position);
	FEs=FEs+1;    
%     fitness_history(i,it)=pop(i).fitness;
%     position_history(i,it,:)=pop(i).Position;
%     Trajectories(i,it)=pop(i).Position(1);
%     it=it+1;


end
% 按照函数值排序
[~, SortOrder]=sort([pop.fitness]);
pop=pop(SortOrder);
Bestsolution=pop(1);                    % 保存一次迭代的最佳值的解决方案
% Convergence_curve=zeros(MaxFEs,1);            % 保留每一次迭代的最佳值
w=1/(sqrt(2*pi)*q*k)*exp(-0.5*(((1:k)-1)/(q*k)).^2);            %计算权重
p_g=w/sum(w);                 % 计算概率

N=SearchAgents_no;


%% 主循环
 while FEs<=Max_iteration   
    % 原k个方案
    s=zeros(k,dim);
    for l=1:k
        s(l,:)=pop(l).Position(1,:);
    end  
    % 新m个方案
    newpop=repmat(empty_individual,m,1);
    for t=1:m
        % 初始化新的解决方案
        newpop(t).Position=zeros(solution);        
        % 生成新的解决方案
        for i=1:dim 
            % 轮盘赌选择高斯函数
            for j=1:k
                if p_g(j)>=rand
                      l=j;
                      break;
                end
            end     
            % 标准差计算
            D=0;
            for r=1:k
                D=D+abs(s(l,i)-s(r,i));
            end
            sigma=ibslo*D/(k-1);           
            % 根据高斯核函数产生高斯随机值
%             newpop(t).Position(i)=normrnd(s(l,i),sigma);  
            newpop(t).Position(i)=s(l,i)+sigma*randn;
        end
        Flag4ub=newpop(t).Position>ub;        %返回大于ub的逻辑值
		Flag4lb=newpop(t).Position<lb;
		newpop(t).Position=(newpop(t).Position.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; 
        % Evaluation
        newpop(t).fitness=fobj(newpop(t).Position); 
        FEs=FEs+1;
        
    end    
    % 将新生成的m个解加入k个原始方案中
    pop=[pop;newpop];
    

    % 将所有解排序
    [~, SortOrder]=sort([pop.fitness]);
    pop=pop(SortOrder);
    % 删除m个解
    pop=pop(1:k);
    
   %% CrissCross
    for l=1:k
        X(l,:)=pop(l).Position(1,:);
        fitness(l)=pop(l).fitness;
    end
    for l=1:k
        pop(l).Position(1,:)=X(l,:);
        pop(l).fitness=fitness(l);
    end  
   

    % 更新最佳值的解决方案
    Bestsolution=pop(1);
    % 保存每次迭代的最佳函数值
    Convergence_curve(it)=Bestsolution.fitness;
    Destination_fitness=Bestsolution.fitness;
    Leader_pos=Bestsolution.Position;
    % 显示每一次的迭代结果信息
%     display(['Iteration ' num2str(it) ': Best fitness = ' num2str(Convergence_curve(it))]);
    it=it+1;
 end
end
 
