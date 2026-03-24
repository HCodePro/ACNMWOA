function [Leader_pos,Convergence_curve]=RCACO(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
solution=[1 dim];             %初始化一个解决方法的大小
%% 初始化相关参数 
FEs=0;
it=1;
k=10;                           % 档案袋大小k
m=SearchAgents_no;              % 新生成个数
q=0.8;     
ibslo=0.9;             
%% 初始化
% 初始化一个个体
empty_individual.Position=[];
empty_individual.fitness=[];
% 初始化种群
pop=repmat(empty_individual,k,1);
for i=1:k
    pop(i).Position=unifrnd(lb,ub,solution);
    Flag4ub=pop(i).Position>ub;        				%返回大于ub的逻辑值
    Flag4lb=pop(i).Position<lb;
    pop(i).Position=(pop(i).Position.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
    pop(i).fitness=fobj(pop(i).Position);
    FEs=FEs+1;
end
% 按照函数值排序
[~, SortOrder]=sort([pop.fitness]);
pop=pop(SortOrder);
w=1/(sqrt(2*pi)*q*k)*exp(-0.5*(((1:k)-1)/(q*k)).^2);            %计算权重
p=w/sum(w);                 % 计算概率
%% 主循环
 while FEs<=MaxFEs   
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
				if p(j)>=rand
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
			newpop(t).Position(i)=s(l,i)+sigma*randn;
		end
		Flag4ub=newpop(t).Position>ub;        %返回大于ub的逻辑值
		Flag4lb=newpop(t).Position<lb;
		newpop(t).Position=(newpop(t).Position.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; 
		
		
		%% 替换机制
		M=newpop(t).Position;
		for h=1:dim
			if(tan(pi*(rand-0.5))<(1-FEs/MaxFEs))  %根据算法剩余运行次数占总运行次数的比值与柯西随机数相比较，使当前位置有一定几率向最优位置靠拢，越后期替换概率越小
				M(h)= pop(1).Position(h);               
			end
		end   
		Fitnessm=fobj(M);               %计算适应度值
		FEs=FEs+1;
		if (Fitnessm<pop(1).fitness)
		   pop(1).fitness = Fitnessm;
		   pop(1).Position =M;
		end 
		%% Evaluation
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
	
	%% CLS机制
    [ pop(1).Position,pop(1).fitness,FEs] = CLS(m,MaxFEs,FEs,lb,ub,pop(1).Position,pop(1).fitness,fobj ); %N：种群大小，Max_iteration：最大评估次数，FEs：当前迭代次数，lb：下界，ub:上届，Destination_position：目标解，Destination_fitness：目标解的函数值，fobj:目标函数 
    %% 
    % 更新最佳值的解决方案
    Bestsolution=pop(1);
    % 保存每次迭代的最佳函数值
    Convergence_curve(it)=Bestsolution.fitness;
    Leader_pos=Bestsolution.Position;
    % 显示每一次的迭代结果信息
%     display(['Iteration ' num2str(it) ': Best fitness = ' num2str(Convergence_curve(it))]);
    it=it+1;
 end
