function [best_pos,Convergence_curve]=MOFOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
tic
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs=0;
it =1;
M = 3;   %子群数||确保各子群数为整数即可
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
X_axis = sx*rand(M,dim);
Y_axis = sy*rand(M,dim);
X = zeros(sizepop/M, dim,M);
Y = zeros(sizepop/M, dim,M);
D = ones(sizepop/M, dim,M);
S = zeros(sizepop/M, dim,M);
Smell = zeros(sizepop/M, 1,M); 
S_temp = zeros(sizepop/M, dim,M);
bestCVaccuarcy = zeros(M,1);

%% 果蝇寻优开始
%利用嗅觉寻找食物
for s = 1:M
for i = 1:sizepop/M
    for j = 1:dim
        % 在全局范围内初始化种群个体
          S(i,j,s) = min(lb) + (max(ub) - min(lb)) * rand();
    end
    
    S_temp(i,:,s) = S(i,:,s);
    %纠正越界个体
    Flag4ub=S(i,:,s)>ub;
    Flag4lb=S(i,:,s)<lb;
    S(i,:,s)=(S(i,:,s).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; 

    FEs = FEs+1;
    Smell(i,1,s) = fobj(S(i,1:dim,s));
end

%根据初始味道浓度寻找初始极值
[bestSmell, bestindex] = min(Smell(:,1,s));
X_axis(s,:) = X(bestindex,:,s);
Y_axis(s,:) = Y(bestindex,:,s);
best_pos = S(bestindex,:,s);
bestCVaccuarcy(s) = bestSmell;
end


%果蝇迭代寻优开始

 while FEs<MaxFEs   
     
    a=2-FEs*((2)/MaxFEs);   %线性递减变量，后期可尝试非线性递减对熟练速度和精度的影响
    for s = 1:M           %离群搜索多种群策略开始【该机制对收敛速度效果不显著，但对混合固定维函数的精度有显著的提升】
    for i = 1:sizepop/M
        
        tempx = ax*rand() - bx;     
        tempy = ay*rand() - by;
        for j = 1:dim
            
            r1 = rand();
            oumiga = min(lb)+(max(ub)-min(lb))*r1;  %添加全局扰动变量，提高子群全局搜索能力
        
            X(i,j,s) = X_axis(s,j) + tempx;
            Y(i,j,s) = Y_axis(s,j) + tempy;
            D(i,j,s) = (X(i,j,s).^2 + Y(i,j,s).^2).^0.5;
            S_temp(i,j,s) = 1./D(i,j,s);
            d = ceil(rand()*dim);      %随机挑选子群中的个体进行搜索半径随迭代数次数递减的全局搜索
            if j == d
                S_temp(i,j,s) = S_temp(i,j,s) + sign(rand() - 0.5) * oumiga*a;
            end

        end
        %出界调整
        if S_temp(i,:,s) > ub
            S_temp(i,:,s) = ub;
        end
        if S_temp(i,:,s) < lb
            S_temp(i,:,s) = lb;
        end
        %前哨果蝇预先探索，判断适应度值是否值得大部队前进。若值得则按照探索方向前进，否则按照原路径前进
        if FEs<MaxFEs
            FEs = FEs+2;
            fitness_st = fobj(S(i,:,s));
            fitness_temp = fobj(S_temp(i,:,s));
            fitness_comp = [fitness_temp,fitness_st];
            [~,p] = min(fitness_comp);
            if p==1
                S(i,:,s) = S_temp(i,:,s);
            else
                S(i,:,s) = S(i,:,s);
            end 
        end
  
        %高斯分布映射以及适应度值选择 【注意：高斯映射机制和前哨搜索机制两者必须同时存在才有显著效果，两者中只存在任意一种均会使效果失效】
        %（高斯映射机制和前哨搜索机制对单峰函数和多峰函数的收敛速度有显著的提升，但对混合固定维函数的效果不太好）
        FEs = FEs+2;
        if FEs<MaxFEs
            Moth_pos_m_gaus=S(i,:,s)*(1+randn(1));
            Moth_fitness_m_gaus=fobj(Moth_pos_m_gaus);
            Moth_fitness_s=fobj(S(i,:,s));
            Moth_fitness_comb=[Moth_fitness_m_gaus,Moth_fitness_s];
            [~,m]=min(Moth_fitness_comb);
            if m==1    
                S(i,:,s)=Moth_pos_m_gaus;
            else
                S(i,:,s)=S(i,:,s); 
            end    
            %越界调整
            Flag4ub=S(i,:,s)>ub;
            Flag4lb=S(i,:,s)<lb;
            S(i,:,s)=(S(i,:,s).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;    
            Smell(i,1,s) = Moth_fitness_s;
        else
            break;
        end


    end  
    
	%根据味道浓度寻找极值
	[bestSmell, bestindex] = min(Smell(:,1,s));
    [~,worstindex] = max(Smell(:,1,s));
    
 	%迭代保留最佳值位置与味道浓度
    %若下一代气味浓度更佳，则按照最佳位置更新气味浓度值和位置，否则当前位置减去最差位置，简称【趋利避害】（但其实此部分优化只对函数5有明显效果，对其他函数目前没观察出显著效果）
	if bestSmell < bestCVaccuarcy
        X_axis(s,:) = X_axis(s,:)+X(bestindex,:,s);
        Y_axis(s,:) = Y_axis(s,:)+Y(bestindex,:,s);
        best_pos=S(bestindex,:,s);
        bestCVaccuarcy(s) = bestSmell;
    else
        X_axis(s,:) = X_axis(s,:)-X(worstindex,:,s);
        Y_axis(s,:) = Y_axis(s,:)-Y(worstindex,:,s);
    end
    
    end
     Convergence_curve(it) = min(bestCVaccuarcy);
     it = it+1;
 end
toc
end

